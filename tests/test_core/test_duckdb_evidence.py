"""M5 acceptance test: evidence cohort ingestion into DuckDB."""

from getRPF.core.release import build_evidence, write_evidence
from getRPF.duckdb.ingest import ingest_cohort


def test_ingest_cohort_batch_loads_evidence_directory(tmp_path):
    evidence_dir = tmp_path / "evidence"
    for sample_id, output in [("s1", 1000), ("s2", 5)]:
        e = build_evidence(
            sample_id=sample_id,
            getrpf_version="0.3.0",
            read_count_input=1000,
            read_count_output=output,
            biological_screen={
                "length_shape": "peaked",
                "insert_mode": 30,
                "adapter_dimer_fraction": 0.0,
                "contamination_screened": False,
                "contamination": {},
                "duplication_umi_aware": 0.1,
                "verdict": "consistent_with_riboseq",
                "reasons": ["peaked"],
            },
        )
        write_evidence(e, evidence_dir / f"{sample_id}.evidence.json")

    db_path = tmp_path / "qc.duckdb"
    ingest_cohort(db_path, evidence_dir)

    import duckdb

    con = duckdb.connect(str(db_path))
    rows = con.execute(
        "SELECT sample_id, release_class FROM evidence ORDER BY sample_id"
    ).fetchall()
    con.close()

    assert rows == [("s1", "clean_no_trim"), ("s2", "exclude_or_hold")]
