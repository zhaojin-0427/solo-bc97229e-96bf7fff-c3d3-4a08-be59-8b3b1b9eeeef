"""Job lifecycle management: persistence, background execution and recovery."""
from __future__ import annotations

import json
import os
import shutil
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .config import Settings, settings as default_settings
from .pipeline import ARTIFACTS, Pipeline, PipelineError

TERMINAL_STATUSES = ("completed", "failed")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class JobManager:
    def __init__(self, settings: Settings = default_settings) -> None:
        self.settings = settings
        self.jobs_dir = settings.jobs_dir
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._state_lock = threading.Lock()
        # In-memory throttle timestamps for progress writes; never persisted.
        self._progress_last_written: dict[str, float] = {}
        self._executor: Optional[ThreadPoolExecutor] = None
        self._executor_lock = threading.Lock()

    # ------------------------------------------------------------------ setup

    def start(self) -> None:
        """Create the worker pool and recover jobs interrupted by a restart."""
        with self._executor_lock:
            if self._executor is None:
                self._executor = ThreadPoolExecutor(
                    max_workers=self.settings.max_workers,
                    thread_name_prefix="fastq-job",
                )
        self.recover_interrupted()

    def shutdown(self) -> None:
        with self._executor_lock:
            if self._executor is not None:
                self._executor.shutdown(wait=False, cancel_futures=True)
                self._executor = None

    def recover_interrupted(self) -> int:
        """Mark every non-terminal job as failed after a process restart."""
        recovered = 0
        if not self.jobs_dir.exists():
            return recovered
        for entry in self.jobs_dir.iterdir():
            if not entry.is_dir():
                continue
            status_path = entry / "status.json"
            if not status_path.exists():
                continue
            try:
                state = self._read_state(entry.name)
            except (OSError, json.JSONDecodeError):
                continue
            if state["status"] not in TERMINAL_STATUSES:
                state["status"] = "failed"
                state["finished_at"] = _utc_now()
                state["error"] = {
                    "message": (
                        "service was interrupted before the job finished; "
                        "the job has been marked failed"
                    ),
                    "side": None,
                    "line_number": None,
                    "record_number": None,
                    "compressed_offset": None,
                }
                self._write_state(entry.name, state)
                self._remove_staging(entry / "work")
                recovered += 1
        return recovered

    # -------------------------------------------------------------- lifecycle

    def create_job(self) -> dict:
        job_id = uuid.uuid4().hex
        job_dir = self.jobs_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=False)
        (job_dir / "inputs").mkdir()
        (job_dir / "work").mkdir()
        state = {
            "job_id": job_id,
            "status": "uploading",
            "created_at": _utc_now(),
            "started_at": None,
            "finished_at": None,
            "error": None,
            "progress": {
                "stage": "uploading",
                "percent": 0,
                "r1_records_parsed": 0,
                "r2_records_parsed": 0,
                "r1_compressed_bytes": 0,
                "r2_compressed_bytes": 0,
                "paired_count": 0,
                "r1_orphan_count": 0,
                "r2_orphan_count": 0,
            },
            "result": None,
            "artifacts": {},
        }
        self._write_state(job_id, state)
        return state

    def submit(self, job_id: str) -> None:
        state = self._read_state(job_id)
        if state["status"] != "uploading":
            raise RuntimeError(f"job {job_id} is not in uploading state")
        state["status"] = "running"
        state["started_at"] = _utc_now()
        self._write_state(job_id, state)
        assert self._executor is not None, "JobManager.start() was not called"
        self._executor.submit(self._run_job, job_id)

    def _run_job(self, job_id: str) -> None:
        job_dir = self._job_dir(job_id)
        try:
            pipeline = Pipeline(
                job_dir=job_dir / "work",
                r1_path=job_dir / "inputs" / "R1.fastq.gz",
                r2_path=job_dir / "inputs" / "R2.fastq.gz",
                batch_size=self.settings.insert_batch_size,
                sqlite_cache_kib=self.settings.sqlite_cache_kib,
                gzip_level=self.settings.gzip_level,
                progress=lambda p: self._update_progress(job_id, p),
            )
            manifest = pipeline.run()
            state = self._read_state(job_id)
            state["status"] = "completed"
            state["finished_at"] = _utc_now()
            state["result"] = {
                "r1_total": manifest["r1_total"],
                "r2_total": manifest["r2_total"],
                "paired_count": manifest["paired_count"],
                "r1_orphan_count": manifest["r1_orphan_count"],
                "r2_orphan_count": manifest["r2_orphan_count"],
            }
            state["artifacts"] = {
                key: filename for key, filename in ARTIFACTS
            }
            state["artifacts"]["manifest"] = "manifest.json"
            state["progress"]["stage"] = "completed"
            state["progress"]["percent"] = 100
            self._write_state(job_id, state)
        except PipelineError as exc:
            self._fail(job_id, exc.to_detail())
            self._remove_staging(job_dir / "work")
        except Exception as exc:  # noqa: BLE001 - surface every failure in state
            self._fail(
                job_id,
                {
                    "message": f"unexpected error: {type(exc).__name__}: {exc}",
                    "side": None,
                    "line_number": None,
                    "record_number": None,
                    "compressed_offset": None,
                },
            )
            self._remove_staging(job_dir / "work")

    def _fail(self, job_id: str, error_detail: dict) -> None:
        state = self._read_state(job_id)
        state["status"] = "failed"
        state["finished_at"] = _utc_now()
        state["error"] = error_detail
        self._write_state(job_id, state)

    def mark_failed_upload(self, job_id: str, message: str) -> None:
        self._fail(
            job_id,
            {
                "message": message,
                "side": None,
                "line_number": None,
                "record_number": None,
                "compressed_offset": None,
            },
        )

    # --------------------------------------------------------------- progress

    def _update_progress(self, job_id: str, progress: dict) -> None:
        # Throttle state writes to roughly one per second; the publishing
        # update is always written through.
        with self._state_lock:
            state_path = self._job_dir(job_id) / "status.json"
            try:
                with open(state_path, "r", encoding="utf-8") as fh:
                    state = json.load(fh)
            except (OSError, json.JSONDecodeError):
                return
            now = datetime.now(timezone.utc).timestamp()
            last = self._progress_last_written.get(job_id)
            if (
                last is not None
                and now - last < 1.0
                and progress["stage"] != "publishing"
            ):
                return
            state["progress"] = progress
            self._atomic_write(state_path, state)
            self._progress_last_written[job_id] = now

    # ----------------------------------------------------------------- queries

    def get_job(self, job_id: str) -> Optional[dict]:
        if not self._is_safe_id(job_id):
            return None
        status_path = self._job_dir(job_id) / "status.json"
        if not status_path.exists():
            return None
        return self._read_state(job_id)

    def list_jobs(self) -> list[dict]:
        if not self.jobs_dir.exists():
            return []
        jobs: list[dict] = []
        for entry in sorted(self.jobs_dir.iterdir(), reverse=True):
            if entry.is_dir() and (entry / "status.json").exists():
                try:
                    jobs.append(self._read_state(entry.name))
                except (OSError, json.JSONDecodeError):
                    continue
        return jobs

    def artifact_path(self, job_id: str, artifact_key: str) -> Optional[Path]:
        state = self.get_job(job_id)
        if state is None or state["status"] != "completed":
            return None
        filename = state.get("artifacts", {}).get(artifact_key)
        if filename is None:
            return None
        path = (self._job_dir(job_id) / "work" / "outputs" / filename).resolve()
        outputs_root = (self._job_dir(job_id) / "work" / "outputs").resolve()
        # Defense in depth against path traversal.
        if outputs_root not in path.parents and path != outputs_root:
            return None
        if not path.is_file():
            return None
        return path

    # ------------------------------------------------------------- spooling I/O

    def spool_upload(
        self, job_id: str, side: str, source, *, content_length: Optional[int]
    ) -> Path:
        """Stream an uploaded gzip into the job's input directory."""
        target = self._job_dir(job_id) / "inputs" / f"{side}.fastq.gz"
        written = 0
        with open(target, "wb") as fh:
            while True:
                chunk = source.read(self.settings.upload_chunk_size)
                if not chunk:
                    break
                fh.write(chunk)
                written += len(chunk)
            fh.flush()
            os.fsync(fh.fileno())
        if content_length is not None and written != content_length:
            raise IOError(
                f"truncated upload for {side}: received {written} "
                f"of {content_length} bytes"
            )
        return target

    # ------------------------------------------------------------- internals

    @staticmethod
    def _is_safe_id(job_id: str) -> bool:
        try:
            uuid.UUID(hex=job_id)
        except (ValueError, AttributeError):
            return False
        return len(job_id) == 32

    def _job_dir(self, job_id: str) -> Path:
        if not self._is_safe_id(job_id):
            raise ValueError("invalid job id")
        return self.jobs_dir / job_id

    def _read_state(self, job_id: str) -> dict:
        with self._state_lock:
            with open(self._job_dir(job_id) / "status.json", "r", encoding="utf-8") as fh:
                return json.load(fh)

    def _write_state(self, job_id: str, state: dict) -> None:
        with self._state_lock:
            self._atomic_write(self._job_dir(job_id) / "status.json", state)

    @staticmethod
    def _atomic_write(path: Path, payload: dict) -> None:
        tmp = path.with_name(path.name + ".tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.rename(tmp, path)

    @staticmethod
    def _remove_staging(work_dir: Path) -> None:
        staging = work_dir / "outputs.staging"
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)
