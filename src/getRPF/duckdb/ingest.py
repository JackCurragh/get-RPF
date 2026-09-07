"""DuckDB ingestion utilities for getRPF reports.

Collects per-sample metrics from:
  - getRPF cleanliness reports (enhanced and basic .rpf_checks)
  - STAR alignment JSON (align-detect)
  - Alignment-based extraction report JSON

Creates/updates a DuckDB database with three normalized tables and a
lightweight view for QC dashboards and gating.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import duckdb

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS samples (
  sample_id TEXT PRIMARY KEY,
  source_path TEXT,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS getrpf_checks (
  sample_id TEXT,
  report_path TEXT,
  is_clean BOOLEAN,
  primary_failure TEXT,
  failure_categories TEXT,
  length_fraction_in_range DOUBLE,
  mode_length INTEGER,
  gc_fraction DOUBLE,
  UNIQUE(sample_id)
);

CREATE TABLE IF NOT EXISTS alignment_stats (
  sample_id TEXT,
  input_reads BIGINT,
  aligned_reads BIGINT,
  alignment_rate DOUBLE,
  mean_5prime_clips DOUBLE,
  mean_3prime_clips DOUBLE,
  recommended_5prime_trim INTEGER,
  recommended_3prime_trim INTEGER,
  consensus_level DOUBLE,
  UNIQUE(sample_id)
);

CREATE TABLE IF NOT EXISTS extraction_summary (
  sample_id TEXT,
  input_reads BIGINT,
  extracted_rpfs BIGINT,
  extraction_rate DOUBLE,
  method TEXT,
  UNIQUE(sample_id)
);

CREATE TABLE IF NOT EXISTS evidence (
  sample_id TEXT,
  release_class TEXT,
  biological_confirmation TEXT,
  retained_fraction DOUBLE,
  read_count_input BIGINT,
  read_count_output BIGINT,
  biological_screen_verdict TEXT,
  evidence_path TEXT,
  evidence_json TEXT,
  UNIQUE(sample_id)
);

CREATE VIEW IF NOT EXISTS qc_overview AS
SELECT s.sample_id,
       g.is_clean,
       g.primary_failure,
       a.alignment_rate,
       a.mean_5prime_clips + a.mean_3prime_clips AS mean_total_clips,
       a.recommended_5prime_trim,
       a.recommended_3prime_trim,
       e.extracted_rpfs,
       e.extraction_rate,
       v.release_class,
       v.biological_confirmation,
       v.biological_screen_verdict
FROM samples s
LEFT JOIN getrpf_checks g USING (sample_id)
LEFT JOIN alignment_stats a USING (sample_id)
LEFT JOIN extraction_summary e USING (sample_id)
LEFT JOIN evidence v USING (sample_id);
"""


def _ensure_db(db_path: Path) -> duckdb.DuckDBPyConnection:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(db_path))
    con.execute(SCHEMA_SQL)
    return con


def ingest_cleanliness(
    con: duckdb.DuckDBPyConnection,
    sample_id: str,
    report_path: Path,
) -> None:
    """Parse the enhanced cleanliness .txt report written by check-cleanliness.

    We look for the Sample Status and any key metrics if present.
    """
    is_clean = None
    primary_failure = None
    failure_categories = None
    length_fraction_in_range = None
    mode_length = None
    gc_fraction = None

    try:
        text = Path(report_path).read_text(errors="replace")
        for line in text.splitlines():
            if line.startswith("Sample Status:"):
                is_clean = "CLEAN" in line
            elif line.startswith("Primary Failure Type:"):
                primary_failure = line.split(":", 1)[1].strip() or None
            elif line.startswith("All Failure Types:"):
                failure_categories = line.split(":", 1)[1].strip() or None
            elif line.strip().startswith("fraction_in_range:"):
                try:
                    length_fraction_in_range = float(line.split(":", 1)[1])
                except Exception:
                    pass
            elif line.strip().startswith("mode_length:"):
                try:
                    mode_length = int(line.split(":", 1)[1])
                except Exception:
                    pass
            elif line.strip().startswith("gc_content:"):
                try:
                    gc_fraction = float(line.split(":", 1)[1])
                except Exception:
                    pass
    except FileNotFoundError:
        pass

    con.execute(
        """
        INSERT INTO getrpf_checks AS g
        (sample_id, report_path, is_clean, primary_failure, failure_categories,
         length_fraction_in_range, mode_length, gc_fraction)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (sample_id) DO UPDATE SET
          report_path=excluded.report_path,
          is_clean=excluded.is_clean,
          primary_failure=excluded.primary_failure,
          failure_categories=excluded.failure_categories,
          length_fraction_in_range=excluded.length_fraction_in_range,
          mode_length=excluded.mode_length,
          gc_fraction=excluded.gc_fraction;
        """,
        [
            sample_id,
            str(report_path),
            is_clean,
            primary_failure,
            failure_categories,
            length_fraction_in_range,
            mode_length,
            gc_fraction,
        ],
    )


def ingest_alignment_json(
    con: duckdb.DuckDBPyConnection,
    sample_id: str,
    align_json: Path,
) -> None:
    data: Dict[str, Any] = {}
    try:
        data = json.loads(Path(align_json).read_text())
    except Exception:
        pass

    stats = data.get("alignment_statistics", {})
    trims = data.get("trim_recommendations", {})

    con.execute(
        """
        INSERT INTO alignment_stats AS a
        (sample_id, input_reads, aligned_reads, alignment_rate,
         mean_5prime_clips, mean_3prime_clips,
         recommended_5prime_trim, recommended_3prime_trim, consensus_level)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (sample_id) DO UPDATE SET
          input_reads=excluded.input_reads,
          aligned_reads=excluded.aligned_reads,
          alignment_rate=excluded.alignment_rate,
          mean_5prime_clips=excluded.mean_5prime_clips,
          mean_3prime_clips=excluded.mean_3prime_clips,
          recommended_5prime_trim=excluded.recommended_5prime_trim,
          recommended_3prime_trim=excluded.recommended_3prime_trim,
          consensus_level=excluded.consensus_level;
        """,
        [
            sample_id,
            stats.get("input_reads"),
            stats.get("aligned_reads"),
            stats.get("alignment_rate"),
            stats.get("mean_5prime_clips"),
            stats.get("mean_3prime_clips"),
            (trims.get("recommended_5prime_trim") or 0),
            (trims.get("recommended_3prime_trim") or 0),
            trims.get("consensus_level"),
        ],
    )


def ingest_extraction_json(
    con: duckdb.DuckDBPyConnection,
    sample_id: str,
    extract_json: Path,
) -> None:
    data: Dict[str, Any] = {}
    try:
        data = json.loads(Path(extract_json).read_text())
    except Exception:
        pass

    summ = data.get("extraction_summary", data)

    con.execute(
        """
        INSERT INTO extraction_summary AS e
        (sample_id, input_reads, extracted_rpfs, extraction_rate, method)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT (sample_id) DO UPDATE SET
          input_reads=excluded.input_reads,
          extracted_rpfs=excluded.extracted_rpfs,
          extraction_rate=excluded.extraction_rate,
          method=excluded.method;
        """,
        [
            sample_id,
            summ.get("input_reads"),
            summ.get("extracted_rpfs"),
            summ.get("extraction_rate"),
            summ.get("method"),
        ],
    )


def ingest_evidence(
    con: duckdb.DuckDBPyConnection,
    sample_id: str,
    evidence_json: Path,
) -> None:
    """Ingest one sample's M5 evidence.json (section 9 schema)."""
    data: Dict[str, Any] = {}
    try:
        data = json.loads(Path(evidence_json).read_text())
    except Exception:
        pass

    biological_screen = data.get("biological_screen", {})

    con.execute(
        """
        INSERT INTO evidence AS v
        (sample_id, release_class, biological_confirmation, retained_fraction,
         read_count_input, read_count_output, biological_screen_verdict,
         evidence_path, evidence_json)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT (sample_id) DO UPDATE SET
          release_class=excluded.release_class,
          biological_confirmation=excluded.biological_confirmation,
          retained_fraction=excluded.retained_fraction,
          read_count_input=excluded.read_count_input,
          read_count_output=excluded.read_count_output,
          biological_screen_verdict=excluded.biological_screen_verdict,
          evidence_path=excluded.evidence_path,
          evidence_json=excluded.evidence_json;
        """,
        [
            sample_id,
            data.get("release_class"),
            data.get("biological_confirmation"),
            data.get("retained_fraction"),
            data.get("read_count_input"),
            data.get("read_count_output"),
            biological_screen.get("verdict"),
            str(evidence_json),
            json.dumps(data),
        ],
    )


def ingest_cohort(
    db_path: Path, evidence_dir: Path, pattern: str = "*.evidence.json"
) -> Path:
    """Batch-ingest a directory of `{sample}.evidence.json` files (M6's
    samplesheet cohort output) into the DuckDB store. Unlike `ingest_all`,
    this does not require one CLI invocation per sample."""
    evidence_dir = Path(evidence_dir)
    con = _ensure_db(Path(db_path))
    for evidence_path in sorted(evidence_dir.glob(pattern)):
        sample_id = (
            evidence_path.name[: -len(pattern) + 1]
            if pattern.startswith("*")
            else evidence_path.stem
        )
        try:
            data = json.loads(evidence_path.read_text())
            sample_id = data.get("sample_id", sample_id)
        except Exception:
            pass
        upsert_sample(con, sample_id, evidence_path)
        ingest_evidence(con, sample_id, evidence_path)
    con.close()
    return Path(db_path)


def upsert_sample(
    con: duckdb.DuckDBPyConnection, sample_id: str, source_path: Path
) -> None:
    con.execute(
        """
        INSERT INTO samples (sample_id, source_path) VALUES (?, ?)
        ON CONFLICT (sample_id) DO UPDATE SET source_path=excluded.source_path;
        """,
        [sample_id, str(source_path)],
    )


def ingest_all(
    db_path: Path,
    sample_id: str,
    source_path: Path,
    cleanliness_report: Optional[Path] = None,
    alignment_json: Optional[Path] = None,
    extraction_json: Optional[Path] = None,
) -> Path:
    con = _ensure_db(Path(db_path))
    upsert_sample(con, sample_id, source_path)
    if cleanliness_report:
        ingest_cleanliness(con, sample_id, cleanliness_report)
    if alignment_json:
        ingest_alignment_json(con, sample_id, alignment_json)
    if extraction_json:
        ingest_extraction_json(con, sample_id, extraction_json)
    con.close()
    return Path(db_path)
