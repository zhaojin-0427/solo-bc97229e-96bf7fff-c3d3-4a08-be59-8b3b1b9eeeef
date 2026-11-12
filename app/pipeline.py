"""Streaming, disk-backed paired-FASTQ processing pipeline.

Design constraints (see README):

* Inputs may be larger than memory. Records are streamed one at a time from
  gzip; no full read or ID collection is ever held in Python memory.
* Matching is performed by primary-key lookups in a per-job SQLite database
  whose page cache is capped (``PRAGMA cache_size``), so the working set stays
  on disk instead of on the heap.
* Original four-line records are stored verbatim (as gzip-decompressed bytes)
  and copied byte-for-byte into the outputs.
* Any gzip corruption, malformed record or duplicate ID fails the whole job;
  partial outputs are written under a staging directory that is only atomically
  renamed into place after every file has been flushed and fsynced.
"""
from __future__ import annotations

import gzip
import hashlib
import os
import re
import sqlite3
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class PipelineError(Exception):
    """Pipeline failure carrying a precise location in the offending input."""

    def __init__(
        self,
        message: str,
        *,
        side: str,
        line_number: int,
        record_number: int,
        compressed_offset: int,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.side = side
        self.line_number = line_number
        self.record_number = record_number
        self.compressed_offset = compressed_offset

    def to_detail(self) -> dict:
        return {
            "message": self.message,
            "side": self.side,
            "line_number": self.line_number,
            "record_number": self.record_number,
            "compressed_offset": self.compressed_offset,
        }


# ---------------------------------------------------------------------------
# FASTQ validation
# ---------------------------------------------------------------------------

# Header: '@', non-empty ID (no whitespace), terminated by /1 or /2.
_HEADER_RE = re.compile(rb"^@(\S+)/(1|2)\r?\n$")
_SEQ_RE = re.compile(rb"^[ACGTN]+\r?\n$")
# Printable ASCII quality, 0x21..0x7e, same length as the sequence; no CR/LF.
_QUAL_RE = re.compile(rb"^[!-~]+\r?\n$")


@dataclass(frozen=True)
class FastqRecord:
    """One four-line record kept verbatim plus its parsed identity."""

    id: bytes
    side: int  # 1 or 2
    data: bytes  # exactly the four decompressed lines, original line endings


@dataclass
class _ReadLocation:
    line_number: int
    record_number: int
    compressed_offset: int


class _GzipFastqReader:
    """Line-by-line gzip reader that tracks position for error reporting."""

    def __init__(self, path: Path, side: str):
        self.side = side
        self._gz = gzip.GzipFile(filename=str(path), mode="rb")
        self._line_number = 0
        self._record_number = 0

    def _compressed_offset(self) -> int:
        # GzipFile.fileobj is the underlying binary file; its position is the
        # number of compressed bytes consumed.
        fileobj = getattr(self._gz, "fileobj", None)
        if fileobj is not None and hasattr(fileobj, "tell"):
            try:
                return int(fileobj.tell())
            except OSError:
                return -1
        return -1

    def _readline(self) -> bytes:
        try:
            line = self._gz.readline()
        except (OSError, EOFError, zlib.error, gzip.BadGzipFile) as exc:
            raise PipelineError(
                f"gzip decompression failed: {type(exc).__name__}: {exc}",
                side=self.side,
                line_number=self._line_number + 1,
                record_number=self._record_number + 1,
                compressed_offset=self._compressed_offset(),
            ) from exc
        self._line_number += 1
        return line

    def read_record(self) -> Optional[tuple[FastqRecord, _ReadLocation]]:
        """Return the next record, or ``None`` at a clean end of file."""
        start_line = self._line_number + 1
        start_offset = self._compressed_offset()
        lines: list[bytes] = []
        for _ in range(4):
            line = self._readline()
            if not line:
                if lines:
                    raise PipelineError(
                        "truncated record: expected 4 lines but reached end of file",
                        side=self.side,
                        line_number=start_line,
                        record_number=self._record_number + 1,
                        compressed_offset=start_offset,
                    )
                return None
            lines.append(line)
        header, seq, plus, qual = lines
        match = _HEADER_RE.match(header)
        if not match:
            raise PipelineError(
                "malformed header: expected '@ID/1' (R1) or '@ID/2' (R2), "
                "ID non-empty and without whitespace",
                side=self.side,
                line_number=start_line,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        rec_id, side_token = match.group(1), int(match.group(2))
        expected_side = 1 if self.side == "R1" else 2
        if side_token != expected_side:
            raise PipelineError(
                f"header declares /{side_token} but the record was submitted as {self.side}",
                side=self.side,
                line_number=start_line,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        if not _SEQ_RE.match(seq):
            raise PipelineError(
                "malformed sequence: must be a non-empty line containing only ACGTN",
                side=self.side,
                line_number=start_line + 1,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        if plus not in (b"+\n", b"+\r\n"):
            raise PipelineError(
                "malformed separator line: it must contain only '+'",
                side=self.side,
                line_number=start_line + 2,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        if not _QUAL_RE.match(qual):
            raise PipelineError(
                "malformed quality: must be printable ASCII without whitespace",
                side=self.side,
                line_number=start_line + 3,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        seq_len = len(seq.rstrip(b"\r\n"))
        qual_len = len(qual.rstrip(b"\r\n"))
        if seq_len != qual_len:
            raise PipelineError(
                f"sequence/quality length mismatch: {seq_len} vs {qual_len}",
                side=self.side,
                line_number=start_line + 3,
                record_number=self._record_number + 1,
                compressed_offset=start_offset,
            )
        self._record_number += 1
        record = FastqRecord(
            id=rec_id,
            side=expected_side,
            data=b"".join(lines),
        )
        location = _ReadLocation(
            line_number=start_line,
            record_number=self._record_number,
            compressed_offset=start_offset,
        )
        return record, location

    def close(self) -> None:
        # GzipFile stops at the first end-of-stream member and silently
        # ignores trailing bytes; treat any trailing non-whitespace bytes as
        # corruption so damaged inputs always fail.
        try:
            tail = self._gz.read()
        except (OSError, EOFError, zlib.error, gzip.BadGzipFile) as exc:
            raise PipelineError(
                f"gzip decompression failed: {type(exc).__name__}: {exc}",
                side=self.side,
                line_number=self._line_number + 1,
                record_number=self._record_number + 1,
                compressed_offset=self._compressed_offset(),
            ) from exc
        if tail:
            raise PipelineError(
                "gzip stream contains trailing bytes after end-of-stream marker",
                side=self.side,
                line_number=self._line_number + 1,
                record_number=self._record_number + 1,
                compressed_offset=self._compressed_offset(),
            )
        try:
            self._gz.close()
        except (OSError, EOFError, zlib.error, gzip.BadGzipFile) as exc:
            raise PipelineError(
                f"gzip decompression failed: {type(exc).__name__}: {exc}",
                side=self.side,
                line_number=self._line_number + 1,
                record_number=self._record_number + 1,
                compressed_offset=self._compressed_offset(),
            ) from exc

    @property
    def compressed_offset(self) -> int:
        return self._compressed_offset()


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

ProgressCallback = Callable[[dict], None]

ARTIFACTS = (
    ("r1_paired", "R1_paired.fastq.gz"),
    ("r2_paired", "R2_paired.fastq.gz"),
    ("r1_orphan", "R1_orphan.fastq.gz"),
    ("r2_orphan", "R2_orphan.fastq.gz"),
)


class Pipeline:
    def __init__(
        self,
        job_dir: Path,
        r1_path: Path,
        r2_path: Path,
        *,
        batch_size: int = 1000,
        sqlite_cache_kib: int = 4096,
        gzip_level: int = 6,
        progress: Optional[ProgressCallback] = None,
    ) -> None:
        self.job_dir = job_dir
        self.r1_path = r1_path
        self.r2_path = r2_path
        self.batch_size = batch_size
        self.sqlite_cache_kib = sqlite_cache_kib
        self.gzip_level = gzip_level
        self._progress_cb = progress or (lambda _p: None)
        self._last_pct = -1

    # -- progress -----------------------------------------------------------

    def _report(
        self,
        stage: str,
        *,
        r1_records: int,
        r2_records: int,
        r1_offset: int,
        r2_offset: int,
        r1_size: int,
        r2_size: int,
        paired: int,
        r1_orphan: int,
        r2_orphan: int,
    ) -> None:
        total_size = r1_size + r2_size
        if total_size > 0:
            pct = int(
                min(99, 100 * (r1_offset + r2_offset) / total_size)
            )
        else:
            pct = 0
        if stage == "publishing":
            pct = 100
        self._last_pct = pct
        self._progress_cb(
            {
                "stage": stage,
                "percent": pct,
                "r1_records_parsed": r1_records,
                "r2_records_parsed": r2_records,
                "r1_compressed_bytes": r1_offset,
                "r2_compressed_bytes": r2_offset,
                "paired_count": paired,
                "r1_orphan_count": r1_orphan,
                "r2_orphan_count": r2_orphan,
            }
        )

    # -- main entry ---------------------------------------------------------

    def run(self) -> dict:
        db_path = self.job_dir / "pipeline.sqlite3"
        staging = self.job_dir / "outputs.staging"
        final = self.job_dir / "outputs"
        if staging.exists():
            # Leftover from a previous failed attempt in the same job dir.
            self._clean_staging(staging)
        else:
            staging.mkdir()
        if final.exists():
            raise RuntimeError("outputs directory already exists; refusing to overwrite")

        r1_size = self.r1_path.stat().st_size
        r2_size = self.r2_path.stat().st_size
        r1_count = r2_count = paired = 0

        conn = sqlite3.connect(str(db_path))
        try:
            self._configure_conn(conn)
            self._create_schema(conn)

            # ---- ingest R1 ----------------------------------------------
            reader = _GzipFastqReader(self.r1_path, "R1")
            try:
                self._report(
                    "ingesting_r1",
                    r1_records=0, r2_records=0,
                    r1_offset=0, r2_offset=0,
                    r1_size=r1_size, r2_size=r2_size,
                    paired=0, r1_orphan=0, r2_orphan=0,
                )
                pending: list[tuple[bytes, bytes, _ReadLocation]] = []
                while True:
                    read = reader.read_record()
                    if read is None:
                        break
                    record, location = read
                    pending.append((record.id, record.data, location))
                    r1_count += 1
                    if len(pending) >= self.batch_size:
                        self._insert_r1_batch(conn, pending)
                        pending.clear()
                        self._report(
                            "ingesting_r1",
                            r1_records=r1_count, r2_records=0,
                            r1_offset=reader.compressed_offset, r2_offset=0,
                            r1_size=r1_size, r2_size=r2_size,
                            paired=0, r1_orphan=0, r2_orphan=0,
                        )
                if pending:
                    self._insert_r1_batch(conn, pending)
                    pending.clear()
                conn.commit()
            finally:
                reader.close()

            # ---- ingest R2 and match ------------------------------------
            reader = _GzipFastqReader(self.r2_path, "R2")
            try:
                self._report(
                    "ingesting_r2",
                    r1_records=r1_count, r2_records=0,
                    r1_offset=r1_size, r2_offset=0,
                    r1_size=r1_size, r2_size=r2_size,
                    paired=0, r1_orphan=0, r2_orphan=0,
                )
                pending: list[tuple[bytes, bytes, _ReadLocation]] = []
                while True:
                    read = reader.read_record()
                    if read is None:
                        break
                    record, location = read
                    pending.append((record.id, record.data, location))
                    r2_count += 1
                    if len(pending) >= self.batch_size:
                        paired += self._insert_r2_batch(conn, pending)
                        pending.clear()
                        self._report(
                            "ingesting_r2",
                            r1_records=r1_count, r2_records=r2_count,
                            r1_offset=r1_size,
                            r2_offset=reader.compressed_offset,
                            r1_size=r1_size, r2_size=r2_size,
                            paired=paired, r1_orphan=0, r2_orphan=0,
                        )
                if pending:
                    paired += self._insert_r2_batch(conn, pending)
                    pending.clear()
                conn.commit()
            finally:
                reader.close()

            # ---- counts --------------------------------------------------
            r1_orphan = conn.execute(
                "SELECT COUNT(*) FROM r1 WHERE matched = 0"
            ).fetchone()[0]
            r2_orphan = conn.execute(
                "SELECT COUNT(*) FROM r2 WHERE r1_id IS NULL"
            ).fetchone()[0]

            # ---- write outputs ------------------------------------------
            self._report(
                "writing_outputs",
                r1_records=r1_count, r2_records=r2_count,
                r1_offset=r1_size, r2_offset=r2_size,
                r1_size=r1_size, r2_size=r2_size,
                paired=paired, r1_orphan=r1_orphan, r2_orphan=r2_orphan,
            )
            checksums = self._write_outputs(conn, staging)
        finally:
            conn.close()

        # ---- atomic publication ------------------------------------------
        self._report(
            "publishing",
            r1_records=r1_count, r2_records=r2_count,
            r1_offset=r1_size, r2_offset=r2_size,
            r1_size=r1_size, r2_size=r2_size,
            paired=paired, r1_orphan=r1_orphan, r2_orphan=r2_orphan,
        )
        manifest = {
            "r1_total": r1_count,
            "r2_total": r2_count,
            "paired_count": paired,
            "r1_orphan_count": r1_orphan,
            "r2_orphan_count": r2_orphan,
            "files": checksums,
        }
        self._write_manifest(staging, manifest)
        os.rename(staging, final)
        # fsync the job directory so the rename itself is durable.
        self._fsync_dir(self.job_dir)
        return manifest

    # -- database helpers ---------------------------------------------------

    def _configure_conn(self, conn: sqlite3.Connection) -> None:
        conn.execute(f"PRAGMA cache_size = -{int(self.sqlite_cache_kib)}")
        conn.execute("PRAGMA temp_store = FILE")
        conn.execute("PRAGMA journal_mode = DELETE")
        conn.execute("PRAGMA synchronous = FULL")

    @staticmethod
    def _create_schema(conn: sqlite3.Connection) -> None:
        # WITHOUT ROWID keyed tables: primary-key lookups (matching and
        # duplicate detection) are B-tree index lookups, no full table scan
        # and no in-memory ID dictionary.
        conn.execute(
            "CREATE TABLE r1 ("
            "id BLOB PRIMARY KEY, record BLOB NOT NULL, "
            "matched INTEGER NOT NULL DEFAULT 0) WITHOUT ROWID"
        )
        conn.execute(
            "CREATE TABLE r2 ("
            "id BLOB PRIMARY KEY, record BLOB NOT NULL, "
            "r1_id BLOB REFERENCES r1(id)) WITHOUT ROWID"
        )

    def _insert_r1_batch(
        self,
        conn: sqlite3.Connection,
        batch: list[tuple[bytes, bytes, _ReadLocation]],
    ) -> None:
        try:
            conn.execute("BEGIN")
            conn.executemany(
                "INSERT INTO r1 (id, record) VALUES (?, ?)",
                [(rec_id, record) for rec_id, record, _loc in batch],
            )
            conn.commit()
        except sqlite3.IntegrityError as exc:
            conn.rollback()
            dup_id, loc = self._locate_duplicate(conn, batch, "r1")
            raise PipelineError(
                f"duplicate R1 id '{dup_id.decode('ascii', 'replace')}': "
                "ids must be unique within each end",
                side="R1",
                line_number=loc.line_number,
                record_number=loc.record_number,
                compressed_offset=loc.compressed_offset,
            ) from exc

    def _insert_r2_batch(
        self,
        conn: sqlite3.Connection,
        batch: list[tuple[bytes, bytes, _ReadLocation]],
    ) -> int:
        paired_now = 0
        try:
            conn.execute("BEGIN")
            for rec_id, record, _loc in batch:
                row = conn.execute(
                    "SELECT 1 FROM r1 WHERE id = ?", (rec_id,)
                ).fetchone()
                if row is not None:
                    conn.execute(
                        "INSERT INTO r2 (id, record, r1_id) VALUES (?, ?, ?)",
                        (rec_id, record, rec_id),
                    )
                    conn.execute(
                        "UPDATE r1 SET matched = 1 WHERE id = ?", (rec_id,)
                    )
                    paired_now += 1
                else:
                    conn.execute(
                        "INSERT INTO r2 (id, record, r1_id) VALUES (?, ?, NULL)",
                        (rec_id, record),
                    )
            conn.commit()
            return paired_now
        except sqlite3.IntegrityError as exc:
            conn.rollback()
            dup_id, loc = self._locate_duplicate(conn, batch, "r2")
            raise PipelineError(
                f"duplicate R2 id '{dup_id.decode('ascii', 'replace')}': "
                "ids must be unique within each end",
                side="R2",
                line_number=loc.line_number,
                record_number=loc.record_number,
                compressed_offset=loc.compressed_offset,
            ) from exc

    @staticmethod
    def _locate_duplicate(
        conn: sqlite3.Connection,
        batch: list[tuple[bytes, bytes, _ReadLocation]],
        table: str,
    ) -> tuple[bytes, _ReadLocation]:
        """Find the offending id and its exact location in the input stream."""
        seen: set[bytes] = set()
        for rec_id, _record, loc in batch:
            if rec_id in seen:
                return rec_id, loc
            seen.add(rec_id)
            row = conn.execute(
                f"SELECT 1 FROM {table} WHERE id = ?", (rec_id,)
            ).fetchone()
            if row is not None:
                return rec_id, loc
        return b"<unknown>", _ReadLocation(
            line_number=-1, record_number=-1, compressed_offset=-1
        )

    # -- output helpers -----------------------------------------------------

    def _write_outputs(
        self, conn: sqlite3.Connection, staging: Path
    ) -> dict[str, dict[str, str]]:
        queries = {
            "r1_paired": "SELECT record FROM r1 WHERE matched = 1 ORDER BY id",
            "r2_paired": (
                "SELECT r2.record FROM r2 "
                "JOIN r1 ON r2.r1_id = r1.id ORDER BY r1.id"
            ),
            "r1_orphan": "SELECT record FROM r1 WHERE matched = 0 ORDER BY id",
            "r2_orphan": "SELECT record FROM r2 WHERE r1_id IS NULL ORDER BY id",
        }
        checksums: dict[str, dict[str, str]] = {}
        for key, filename in ARTIFACTS:
            target = staging / filename
            # Checksum of the record bytes (content identity, independent of
            # gzip framing/compression level).
            content_digest = hashlib.sha256()
            # Checksum of the actual gzip file that gets downloaded.
            file_digest = hashlib.sha256()
            # Cursor iteration is streaming; one record blob is materialized
            # at a time.
            cursor = conn.execute(queries[key])
            with gzip.GzipFile(
                filename=str(target), mode="wb", compresslevel=self.gzip_level
            ) as out:
                for (record,) in cursor:
                    out.write(record)
                    content_digest.update(record)
                out.flush()
                os.fsync(out.fileno())
            with open(target, "rb") as fh:
                while True:
                    chunk = fh.read(1024 * 1024)
                    if not chunk:
                        break
                    file_digest.update(chunk)
            checksums[filename] = {
                "sha256_records": content_digest.hexdigest(),
                "sha256_gzip_file": file_digest.hexdigest(),
            }
        return checksums

    @staticmethod
    def _write_manifest(staging: Path, manifest: dict) -> None:
        import json

        manifest_path = staging / "manifest.json"
        tmp_path = staging / "manifest.json.tmp"
        with open(tmp_path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.rename(tmp_path, manifest_path)
        Pipeline._fsync_dir(staging)

    @staticmethod
    def _clean_staging(staging: Path) -> None:
        # Leftover from a previous failed attempt in the same job directory.
        for child in staging.iterdir():
            child.unlink()
        staging.rmdir()
        staging.mkdir()

    @staticmethod
    def _fsync_dir(path: Path) -> None:
        fd = os.open(str(path), os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
