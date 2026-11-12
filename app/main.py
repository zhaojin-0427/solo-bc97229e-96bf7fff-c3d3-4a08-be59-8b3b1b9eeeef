"""HTTP API for paired-end FASTQ ordering jobs."""
from __future__ import annotations

from typing import Optional

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse

from .config import Settings, settings as default_settings
from .manager import JobManager

ARTIFACT_MEDIA_TYPES = {
    "manifest": "application/json",
}
ARTIFACT_DOWNLOAD_NAMES = {
    key: filename
    for key, filename in (
        ("r1_paired", "R1_paired.fastq.gz"),
        ("r2_paired", "R2_paired.fastq.gz"),
        ("r1_orphan", "R1_orphan.fastq.gz"),
        ("r2_orphan", "R2_orphan.fastq.gz"),
        ("manifest", "manifest.json"),
    )
}


def create_app(settings: Optional[Settings] = None) -> FastAPI:
    settings = settings or default_settings
    manager = JobManager(settings)

    app = FastAPI(
        title="Paired-end FASTQ Ordering Service",
        version="1.0.0",
        description=(
            "Upload a pair of gzip-compressed paired-end FASTQ files. "
            "The service streams both ends without loading them into memory, "
            "orders them by ID into matching paired R1/R2 files and orphan "
            "files, and publishes all artifacts atomically. Use the "
            "'Try it out' buttons below to walk through upload → status → "
            "download."
        ),
        lifespan=_lifespan(manager),
    )
    app.state.job_manager = manager
    app.state.settings = settings

    @app.get("/", tags=["meta"])
    def root() -> dict:
        return {
            "service": "paired-fastq-ordering",
            "docs": "/docs",
            "openapi": "/openapi.json",
            "endpoints": [
                "POST /jobs",
                "GET  /jobs",
                "GET  /jobs/{job_id}",
                "GET  /jobs/{job_id}/download/{artifact}",
            ],
        }

    @app.get("/health", tags=["meta"])
    def health() -> dict:
        return {"status": "ok"}

    @app.post(
        "/jobs",
        status_code=201,
        tags=["jobs"],
        summary="Upload an R1/R2 gzip FASTQ pair and start processing",
    )
    async def create_job(
        r1: UploadFile = File(
            ...,
            description="R1 end, gzip FASTQ, headers '@ID/1', four lines per record",
        ),
        r2: UploadFile = File(
            ...,
            description="R2 end, gzip FASTQ, headers '@ID/2', four lines per record",
        ),
    ) -> dict:
        state = manager.create_job()
        job_id = state["job_id"]
        try:
            for upload, side in ((r1, "R1"), (r2, "R2")):
                content_length = _content_length(upload)
                # Blocking disk I/O runs in a worker thread so the event
                # loop stays free for status queries from other clients.
                await run_in_threadpool(
                    manager.spool_upload,
                    job_id,
                    side,
                    upload.file,
                    content_length=content_length,
                )
                await upload.close()
        except Exception as exc:  # noqa: BLE001
            manager.mark_failed_upload(job_id, f"upload failed: {exc}")
            raise HTTPException(
                status_code=400,
                detail={
                    "job_id": job_id,
                    "message": f"upload failed: {exc}",
                },
            )
        manager.submit(job_id)
        return manager.get_job(job_id)

    @app.get(
        "/jobs",
        tags=["jobs"],
        summary="List all jobs (newest first)",
    )
    def list_jobs() -> dict:
        return {"jobs": manager.list_jobs()}

    @app.get(
        "/jobs/{job_id}",
        tags=["jobs"],
        summary="Get job progress, counts and (on failure) error location",
    )
    def get_job(job_id: str) -> dict:
        state = manager.get_job(job_id)
        if state is None:
            raise HTTPException(status_code=404, detail="job not found")
        return state

    @app.get(
        "/jobs/{job_id}/download/{artifact}",
        tags=["jobs"],
        summary=(
            "Download a result artifact: r1_paired, r2_paired, r1_orphan, "
            "r2_orphan or manifest. Only completed jobs are downloadable."
        ),
    )
    def download_artifact(job_id: str, artifact: str) -> FileResponse:
        if artifact not in ARTIFACT_DOWNLOAD_NAMES:
            raise HTTPException(
                status_code=404,
                detail=(
                    f"unknown artifact '{artifact}'; choose one of: "
                    + ", ".join(sorted(ARTIFACT_DOWNLOAD_NAMES))
                ),
            )
        state = manager.get_job(job_id)
        if state is None:
            raise HTTPException(status_code=404, detail="job not found")
        if state["status"] != "completed":
            raise HTTPException(
                status_code=409,
                detail=(
                    f"job is '{state['status']}'; artifacts are only "
                    "published after successful completion"
                ),
            )
        path = manager.artifact_path(job_id, artifact)
        if path is None:
            raise HTTPException(status_code=404, detail="artifact not found")
        media_type = ARTIFACT_MEDIA_TYPES.get(
            artifact, "application/gzip"
        )
        return FileResponse(
            path,
            media_type=media_type,
            filename=ARTIFACT_DOWNLOAD_NAMES[artifact],
        )

    return app


def _lifespan(manager: JobManager):
    from contextlib import asynccontextmanager

    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        manager.start()
        try:
            yield
        finally:
            manager.shutdown()

    return lifespan


def _content_length(upload: UploadFile) -> Optional[int]:
    headers = getattr(upload, "headers", None)
    if headers is None:
        return None
    raw = headers.get("content-length")
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


# Module-level app used by uvicorn (`app.main:app`) and Docker. It is created
# lazily so that importing this module (e.g. in tests with custom settings)
# never touches the default data directory.
_app: Optional[FastAPI] = None


def __getattr__(name: str):
    if name == "app":
        global _app
        if _app is None:
            _app = create_app()
        return _app
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
