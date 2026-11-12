"""Shared pytest fixtures and FASTQ test-data builders."""
from __future__ import annotations

import gzip
import io
import os
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.config import Settings
from app.main import create_app


def make_fastq_gz(records: list[tuple[str, str, str]], side: int) -> bytes:
    """Build gzip bytes for records given as (id, sequence, quality)."""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        for rec_id, seq, qual in records:
            gz.write(f"@{rec_id}/{side}\n".encode())
            gz.write(f"{seq}\n".encode())
            gz.write(b"+\n")
            gz.write(f"{qual}\n".encode())
    return buf.getvalue()


def gunzip_records(data: bytes) -> list[bytes]:
    """Split a gzip FASTQ blob into its raw four-line records."""
    text = gzip.decompress(data)
    lines = text.splitlines(keepends=True)
    assert len(lines) % 4 == 0
    return [b"".join(lines[i : i + 4]) for i in range(0, len(lines), 4)]


@pytest.fixture()
def jobs_dir(tmp_path: Path) -> Path:
    path = tmp_path / "jobs"
    path.mkdir()
    return path


@pytest.fixture()
def client(jobs_dir: Path) -> TestClient:
    settings = Settings(jobs_dir=jobs_dir, max_workers=2)
    app = create_app(settings)
    with TestClient(app) as test_client:
        yield test_client


@pytest.fixture()
def env_jobs_dir(monkeypatch, tmp_path: Path) -> Path:
    path = tmp_path / "jobs"
    path.mkdir()
    monkeypatch.setenv("FASTQ_JOBS_DIR", str(path))
    return path


def upload_pair(
    client: TestClient, r1: bytes, r2: bytes
) -> str:
    response = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", r1, "application/gzip"),
            "r2": ("R2.fastq.gz", r2, "application/gzip"),
        },
    )
    assert response.status_code == 201, response.text
    return response.json()["job_id"]


def wait_terminal(client: TestClient, job_id: str) -> dict:
    # TestClient runs the thread pool inside the same process; the background
    # future normally completes quickly. Poll with a generous cap.
    for _ in range(600):
        state = client.get(f"/jobs/{job_id}").json()
        if state["status"] in ("completed", "failed"):
            return state
    raise AssertionError("job did not reach a terminal state")
