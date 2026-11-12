"""End-to-end tests for the paired FASTQ ordering API and pipeline."""
from __future__ import annotations

import gzip
import io
import json
import os
import resource
import threading

from app.config import Settings
from app.main import create_app
from app.manager import JobManager
from app.pipeline import Pipeline
from fastapi.testclient import TestClient

from conftest import gunzip_records, make_fastq_gz, upload_pair, wait_terminal


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_paired_and_orphan_outputs(client: TestClient) -> None:
    # R1: ids a, b, c (c has no R2 mate); R2: ids a, b, d (d has no R1 mate).
    r1 = make_fastq_gz(
        [
            ("idA", "ACGT", "IIII"),
            ("idB", "ACGTN", "IIIII"),
            ("idC", "NNNN", "!!!!"),
        ],
        side=1,
    )
    r2 = make_fastq_gz(
        [
            ("idB", "ACGT", "IIII"),
            ("idA", "TTTT", "####"),
            ("idD", "ACGT", "IIII"),
        ],
        side=2,
    )
    job_id = upload_pair(client, r1, r2)
    state = wait_terminal(client, job_id)

    assert state["status"] == "completed", state
    assert state["result"] == {
        "r1_total": 3,
        "r2_total": 3,
        "paired_count": 2,
        "r1_orphan_count": 1,
        "r2_orphan_count": 1,
    }
    assert state["progress"]["percent"] == 100

    r1_paired = client.get(f"/jobs/{job_id}/download/r1_paired").content
    r2_paired = client.get(f"/jobs/{job_id}/download/r2_paired").content
    r1_orphan = client.get(f"/jobs/{job_id}/download/r1_orphan").content
    r2_orphan = client.get(f"/jobs/{job_id}/download/r2_orphan").content

    r1_paired_records = gunzip_records(r1_paired)
    r2_paired_records = gunzip_records(r2_paired)
    # Same order by ID on both paired ends, starting with idA.
    assert [r.split(b"\n", 1)[0] for r in r1_paired_records] == [
        b"@idA/1",
        b"@idB/1",
    ]
    assert [r.split(b"\n", 1)[0] for r in r2_paired_records] == [
        b"@idA/2",
        b"@idB/2",
    ]
    # Original four-line records are preserved verbatim.
    assert r1_paired_records[0] == b"@idA/1\nACGT\n+\nIIII\n"
    assert r2_paired_records[1] == b"@idB/2\nACGT\n+\nIIII\n"
    assert [r.split(b"\n", 1)[0] for r in gunzip_records(r1_orphan)] == [
        b"@idC/1"
    ]
    assert [r.split(b"\n", 1)[0] for r in gunzip_records(r2_orphan)] == [
        b"@idD/2"
    ]

    manifest = client.get(f"/jobs/{job_id}/download/manifest").json()
    assert manifest["paired_count"] == 2
    assert set(manifest["files"]) == {
        "R1_paired.fastq.gz",
        "R2_paired.fastq.gz",
        "R1_orphan.fastq.gz",
        "R2_orphan.fastq.gz",
    }


def test_empty_inputs_produce_empty_outputs(client: TestClient) -> None:
    empty = make_fastq_gz([], side=1)
    job_id = upload_pair(client, empty, make_fastq_gz([], side=2))
    state = wait_terminal(client, job_id)
    assert state["status"] == "completed"
    assert state["result"]["paired_count"] == 0
    for artifact in ("r1_paired", "r2_paired", "r1_orphan", "r2_orphan"):
        body = client.get(f"/jobs/{job_id}/download/{artifact}").content
        assert gzip.decompress(body) == b""


# ---------------------------------------------------------------------------
# Validation failures: the whole job fails and no artifacts are downloadable
# ---------------------------------------------------------------------------


def test_duplicate_id_fails_with_location(client: TestClient) -> None:
    r1 = make_fastq_gz(
        [("idA", "ACGT", "IIII"), ("idA", "ACGT", "IIII")], side=1
    )
    r2 = make_fastq_gz([("idA", "ACGT", "IIII")], side=2)
    job_id = upload_pair(client, r1, r2)
    state = wait_terminal(client, job_id)

    assert state["status"] == "failed"
    error = state["error"]
    assert "duplicate R1 id 'idA'" in error["message"]
    assert error["side"] == "R1"
    assert error["record_number"] == 2
    assert error["line_number"] == 5  # second record starts on line 5
    assert error["compressed_offset"] >= 0

    # No partial artifact may be downloaded.
    response = client.get(f"/jobs/{job_id}/download/r1_paired")
    assert response.status_code == 409
    # Staging directory must have been removed; only completed jobs have
    # an outputs directory.
    from app.manager import JobManager as _JM  # noqa: F401

    manager = client.app.state.job_manager
    work_dir = manager.jobs_dir / job_id / "work"
    assert not (work_dir / "outputs.staging").exists()
    assert not (work_dir / "outputs").exists()


def test_duplicate_r2_id_fails(client: TestClient) -> None:
    r1 = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
    r2 = make_fastq_gz(
        [
            ("idA", "ACGT", "IIII"),
            ("idZ", "ACGT", "IIII"),
            ("idZ", "ACGT", "IIII"),
        ],
        side=2,
    )
    job_id = upload_pair(client, r1, r2)
    state = wait_terminal(client, job_id)
    assert state["status"] == "failed"
    assert state["error"]["side"] == "R2"
    assert state["error"]["record_number"] == 3


def test_malformed_records_report_exact_line(client: TestClient) -> None:
    # Bad sequence character on line 2 of the only record.
    bad_seq = b"@idA/1\nACGTX\n+\nIIII\n"
    r1 = gzip.compress(bad_seq)
    r2 = make_fastq_gz([("idA", "ACGT", "IIII")], side=2)
    state = wait_terminal(client, upload_pair(client, r1, r2))
    assert state["status"] == "failed"
    assert state["error"]["line_number"] == 2
    assert "sequence" in state["error"]["message"]

    # Quality/sequence length mismatch.
    bad_len = gzip.compress(b"@idA/1\nACGT\n+\nIII\n")
    state = wait_terminal(client, upload_pair(client, bad_len, r2))
    assert state["status"] == "failed"
    assert state["error"]["line_number"] == 4
    assert "length mismatch" in state["error"]["message"]

    # Wrong side token in header.
    bad_side = gzip.compress(b"@idA/2\nACGT\n+\nIIII\n")
    state = wait_terminal(client, upload_pair(client, bad_side, r2))
    assert state["status"] == "failed"
    assert state["error"]["line_number"] == 1

    # Truncated record (only three lines).
    truncated = gzip.compress(b"@idA/1\nACGT\n+\n")
    state = wait_terminal(client, upload_pair(client, truncated, r2))
    assert state["status"] == "failed"
    assert "truncated" in state["error"]["message"]


def test_corrupt_gzip_fails(client: TestClient) -> None:
    good = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
    corrupted = good[: len(good) // 2] + b"\x00\x00\x00garbage"
    state = wait_terminal(
        client,
        upload_pair(
            client,
            corrupted,
            make_fastq_gz([("idA", "ACGT", "IIII")], side=2),
        ),
    )
    assert state["status"] == "failed"
    assert "gzip" in state["error"]["message"].lower()
    assert state["error"]["side"] == "R1"


def test_trailing_bytes_after_gzip_stream_fail(client: TestClient) -> None:
    good = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
    with_trailing = good + b"not-part-of-the-gzip-stream"
    state = wait_terminal(
        client,
        upload_pair(
            client,
            with_trailing,
            make_fastq_gz([("idA", "ACGT", "IIII")], side=2),
        ),
    )
    assert state["status"] == "failed"
    assert (
        "trailing bytes" in state["error"]["message"]
        or "gzip" in state["error"]["message"].lower()
    )


def test_wrong_end_header_rejected(client: TestClient) -> None:
    # A /1 header uploaded in the R2 slot.
    r1 = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
    wrong_r2 = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
    state = wait_terminal(client, upload_pair(client, r1, wrong_r2))
    assert state["status"] == "failed"
    assert "submitted as R2" in state["error"]["message"]


# ---------------------------------------------------------------------------
# Restart recovery
# ---------------------------------------------------------------------------


def test_interrupted_jobs_fail_on_restart(jobs_dir) -> None:
    settings = Settings(jobs_dir=jobs_dir)
    manager = JobManager(settings)
    state = manager.create_job()
    manager._write_state(
        state["job_id"],
        {**state, "status": "running", "started_at": "2026-01-01T00:00:00+00:00"},
    )
    # Simulate a half-published job left behind by the crash.
    stale_staging = jobs_dir / state["job_id"] / "work" / "outputs.staging"
    stale_staging.mkdir(parents=True)
    (stale_staging / "partial.fastq.gz").write_bytes(b"partial")
    # Simulate process restart: a brand-new manager over the same directory.
    new_manager = JobManager(settings)
    recovered = new_manager.recover_interrupted()
    assert recovered == 1
    recovered_state = new_manager.get_job(state["job_id"])
    assert recovered_state["status"] == "failed"
    assert "interrupted" in recovered_state["error"]["message"]
    assert recovered_state["finished_at"] is not None
    assert not stale_staging.exists()
    assert not (jobs_dir / state["job_id"] / "work" / "outputs").exists()


def test_completed_jobs_survive_restart(jobs_dir) -> None:
    settings = Settings(jobs_dir=jobs_dir)
    manager = JobManager(settings)
    with TestClient(create_app(settings)) as client:
        r1 = make_fastq_gz([("idA", "ACGT", "IIII")], side=1)
        r2 = make_fastq_gz([("idA", "ACGT", "IIII")], side=2)
        job_id = upload_pair(client, r1, r2)
        state = wait_terminal(client, job_id)
        assert state["status"] == "completed"

    # Restart: completed state and artifacts remain downloadable.
    with TestClient(create_app(settings)) as client2:
        again = client2.get(f"/jobs/{job_id}").json()
        assert again["status"] == "completed"
        body = client2.get(f"/jobs/{job_id}/download/r1_paired").content
        assert gunzip_records(body)[0].startswith(b"@idA/1")


# ---------------------------------------------------------------------------
# Concurrency isolation
# ---------------------------------------------------------------------------


def test_concurrent_jobs_use_isolated_directories(client: TestClient) -> None:
    payloads = []
    for suffix in ("X", "Y", "Z"):
        r1 = make_fastq_gz(
            [
                (f"shared{suffix}", "ACGT", "IIII"),
                (f"only1{suffix}", "ACGT", "IIII"),
            ],
            side=1,
        )
        r2 = make_fastq_gz([(f"shared{suffix}", "ACGT", "IIII")], side=2)
        payloads.append((r1, r2))

    job_ids: list[str] = []
    errors: list[Exception] = []

    def worker(payload: tuple[bytes, bytes]) -> None:
        try:
            job_ids.append(upload_pair(client, *payload))
        except Exception as exc:  # pragma: no cover - debugging aid
            errors.append(exc)

    threads = [threading.Thread(target=worker, args=(p,)) for p in payloads]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors

    states = [wait_terminal(client, jid) for jid in job_ids]
    assert all(s["status"] == "completed" for s in states), states
    for state in states:
        assert state["result"]["paired_count"] == 1
        assert state["result"]["r1_orphan_count"] == 1
        assert state["result"]["r2_orphan_count"] == 0

    listing = client.get("/jobs").json()["jobs"]
    assert {j["job_id"] for j in listing} == set(job_ids)
    # Each job lives in its own directory; no shared files.
    jobs_root = client.app.state.job_manager.jobs_dir
    job_dirs = [p for p in os.listdir(jobs_root) if os.path.isdir(os.path.join(jobs_root, p))]
    assert sorted(job_dirs) == sorted(job_ids)


# ---------------------------------------------------------------------------
# Bounded-memory guarantee: large streamed job stays within a small RSS cap
# ---------------------------------------------------------------------------


def test_pipeline_memory_is_bounded(tmp_path) -> None:
    # 200k records (~ tens of MB compressed); RSS must stay well under the
    # size of the input, proving no full read/ID set is held in memory.
    n = 200_000
    r1_path = tmp_path / "R1.fastq.gz"
    r2_path = tmp_path / "R2.fastq.gz"
    work = tmp_path / "work"
    work.mkdir()

    with gzip.open(r1_path, "wb", compresslevel=6) as fh:
        for i in range(n):
            fh.write(f"@id{i:08d}/1\nACGTACGT\n+\nIIIIIIII\n".encode())
    # R2 in reverse order so matching cannot rely on input ordering; every
    # second R2 id has no R1 mate -> exercises orphans at scale too.
    with gzip.open(r2_path, "wb", compresslevel=6) as fh:
        for i in range(n - 1, -1, -1):
            if i % 2 == 0:
                fh.write(f"@id{i:08d}/2\nACGTACGT\n+\nIIIIIIII\n".encode())

    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    # 256 MiB address-space ceiling; SQLite cache itself is capped at 4 MiB.
    cap = 256 * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    try:
        pipeline = Pipeline(
            job_dir=work,
            r1_path=r1_path,
            r2_path=r2_path,
            batch_size=1000,
            sqlite_cache_kib=4096,
        )
        manifest = pipeline.run()
    finally:
        resource.setrlimit(resource.RLIMIT_AS, (soft, hard))

    assert manifest["r1_total"] == n
    assert manifest["r2_total"] == n // 2
    assert manifest["paired_count"] == n // 2
    assert manifest["r1_orphan_count"] == n // 2
    assert manifest["r2_orphan_count"] == 0
    outputs = work / "outputs"
    assert outputs.is_dir()
    assert not (work / "outputs.staging").exists()
