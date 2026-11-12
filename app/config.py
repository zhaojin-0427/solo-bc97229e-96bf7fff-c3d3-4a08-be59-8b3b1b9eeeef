"""Runtime configuration, all values overridable through environment variables."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def _int_env(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None or raw.strip() == "":
        return default
    return int(raw)


@dataclass(frozen=True)
class Settings:
    # Root directory holding one sub-directory per job.
    jobs_dir: Path = Path(os.environ.get("FASTQ_JOBS_DIR", "/data/jobs"))
    # Number of jobs processed concurrently in the background thread pool.
    max_workers: int = _int_env("FASTQ_MAX_WORKERS", 2)
    # Upload spool chunk size (bytes) and pipeline batch/commit size (records).
    upload_chunk_size: int = _int_env("FASTQ_UPLOAD_CHUNK", 1024 * 1024)
    insert_batch_size: int = _int_env("FASTQ_INSERT_BATCH", 1000)
    # SQLite page-cache limit in KiB (negative value => KiB). Keeps the DB
    # working set off the Python heap so no full read/ID set lives in memory.
    sqlite_cache_kib: int = _int_env("FASTQ_SQLITE_CACHE_KIB", 4096)
    # gzip output level.
    gzip_level: int = _int_env("FASTQ_GZIP_LEVEL", 6)


settings = Settings()
