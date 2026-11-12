# Paired-end FASTQ Ordering Service

A FastAPI service that takes a pair of gzip-compressed paired-end FASTQ files
(R1/R2), and produces ID-ordered, mutually aligned paired files plus orphan
files for each end — designed for inputs **larger than memory**.

## What it does

For each submitted pair it writes four gzip FASTQ artifacts and a manifest:

| Artifact key      | File                  | Contents                                             |
|-------------------|-----------------------|------------------------------------------------------|
| `r1_paired`       | `R1_paired.fastq.gz`  | R1 records whose ID also exists in R2, ID-sorted     |
| `r2_paired`       | `R2_paired.fastq.gz`  | matching R2 records in **the same ID order**         |
| `r1_orphan`       | `R1_orphan.fastq.gz`  | R1-only IDs, ID-sorted                               |
| `r2_orphan`       | `R2_orphan.fastq.gz`  | R2-only IDs, ID-sorted                               |
| `manifest`        | `manifest.json`       | counts and SHA-256 (record content and gzip file) of every FASTQ artifact |

Records are copied **verbatim** (original four lines, original line endings).
Validation rules:

* exactly four lines per record;
* header `@ID/1` (R1) or `@ID/2` (R2), non-empty ID without whitespace;
* third line is exactly `+`;
* sequence contains only `ACGTN` and is non-empty;
* quality is printable ASCII (0x21–0x7E), same length as the sequence;
* IDs are unique within each end (duplicate fails the whole job);
* gzip must decompress cleanly (truncated/corrupt streams fail the whole job).

## Streaming / bounded-memory design

No full read set or ID set is ever held in Python memory:

1. **Parsing** — `gzip.GzipFile.readline()` feeds a four-line state machine;
   only one record exists in memory at a time. Uploads are spooled in fixed
   1 MiB chunks.
2. **Matching** — every job gets its own SQLite database (`WITHOUT ROWID`
   tables keyed by the record ID). Matching is a primary-key B-tree lookup and
   duplicate detection is the primary-key constraint. The SQLite page cache is
   capped with `PRAGMA cache_size = -4096` (~4 MiB) and temp storage is forced
   to disk, so arbitrarily large jobs stay off the heap.
3. **Writing** — outputs are streamed from SQLite cursors one record at a time
   into gzip writers; nothing is buffered wholesale.
4. **Atomic publication** — artifacts are first written to
   `work/outputs.staging/`, fsynced, then the whole directory is renamed to
   `work/outputs/` in one atomic operation and the parent directory is
   fsynced. Downloads are only served from `completed` jobs, so partial
   artifacts can never be fetched.

A 200 000-record job is exercised under a 256 MiB address-space cap in
`tests/test_pipeline.py::test_pipeline_memory_is_bounded`.

## Failure, restart and concurrency semantics

* Any gzip corruption, malformed record or duplicate ID fails the entire job;
  the error object reports the end (`R1`/`R2`), 1-based line number, record
  number and compressed-stream byte offset.
* Staging output is deleted on failure; no `outputs/` directory is ever
  published for a failed job, and downloads return HTTP 409.
* Job state is written atomically (`status.json.tmp` → rename, fsynced). On
  startup every job not already in a terminal state (`completed`/`failed`) is
  marked failed with an "interrupted" reason. Completed jobs and their
  artifacts are preserved and remain downloadable.
* Each job gets a UUID directory (`jobs/<uuid>/{inputs,work,status.json}`);
  concurrent jobs never share files, and job IDs are validated before any
  path access.

## Run with Docker Compose

```bash
docker compose up --build
# API:        http://localhost:8000
# Swagger UI: http://localhost:8000/docs  (interactive: upload, query, download)
# health:     http://localhost:8000/health
```

Job data lives in the named volume `fastq-jobs-data` (`/data/jobs` in the
container), so restarts preserve completed work and allow interrupted-job
recovery. Tunables (environment variables):

| Variable                  | Default        | Meaning                              |
|---------------------------|----------------|--------------------------------------|
| `FASTQ_JOBS_DIR`          | `/data/jobs`   | per-job storage root                 |
| `FASTQ_MAX_WORKERS`       | `2`            | concurrent background jobs           |
| `FASTQ_INSERT_BATCH`      | `1000`         | records per SQLite transaction       |
| `FASTQ_SQLITE_CACHE_KIB`  | `4096`         | SQLite page-cache cap (KiB)          |
| `FASTQ_GZIP_LEVEL`        | `6`            | gzip compression level for outputs   |
| `FASTQ_UPLOAD_CHUNK`      | `1048576`      | upload spool chunk size (bytes)      |

## Interactive demo through the API docs

1. Open `http://localhost:8000/docs`.
2. **POST /jobs** → *Try it out* → choose `R1.fastq.gz` and `R2.fastq.gz`
   for the `r1` / `r2` form fields → *Execute*. Copy `job_id` from the
   response.
3. **GET /jobs/{job_id}** → paste the id → *Execute*: poll until
   `"status": "completed"` and inspect `result` (paired/orphan counts) and
   `progress`.
4. **GET /jobs/{job_id}/download/{artifact}** → set `artifact` to
   `r1_paired`, `r2_paired`, `r1_orphan`, `r2_orphan` or `manifest` →
   *Execute* → *Download file*.

Equivalent curl flow:

```bash
curl -X POST http://localhost:8000/jobs \
  -F "r1=@R1.fastq.gz;type=application/gzip" \
  -F "r2=@R2.fastq.gz;type=application/gzip"

curl http://localhost:8000/jobs/<job_id>

curl -OJ http://localhost:8000/jobs/<job_id>/download/r1_paired
curl -OJ http://localhost:8000/jobs/<job_id>/download/manifest
```

A failure response looks like:

```json
{
  "status": "failed",
  "error": {
    "message": "duplicate R2 id 'dup': ids must be unique within each end",
    "side": "R2",
    "line_number": 5,
    "record_number": 2,
    "compressed_offset": 61
  }
}
```

## Local development / tests

```bash
pip install -r requirements.txt
pip install pytest httpx
pytest
uvicorn app.main:app --reload
```

## On-disk layout

```
jobs/
└── <uuid>/
    ├── status.json                 # atomically persisted job state
    ├── inputs/
    │   ├── R1.fastq.gz             # spooled uploads
    │   └── R2.fastq.gz
    └── work/
        ├── pipeline.sqlite3        # disk-backed matching tables
        ├── outputs.staging/        # present only while publishing
        └── outputs/                # appears atomically on success
            ├── R1_paired.fastq.gz
            ├── R2_paired.fastq.gz
            ├── R1_orphan.fastq.gz
            ├── R2_orphan.fastq.gz
            └── manifest.json
```
