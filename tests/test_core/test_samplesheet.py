"""M6 acceptance tests: samplesheet driver + audit/override modes.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M6.
"""

from getRPF.core.samplesheet import run_cohort

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def _write_fastq(path, sequences):
    lines = []
    for i, seq in enumerate(sequences):
        lines.append(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    path.write_text("".join(lines))


def _artifact_sequences():
    """20000 reads of length 30 with a 99% 3' 'A' artifact, plus small
    clean control classes -- same shape as the M3 fixture."""
    sequences = []
    for i in range(20000):
        body = _cycled_seq(i, 29)
        tail = "A" if i < 19800 else BASES[i % 4]
        sequences.append(body + tail)
    for length in (28, 29, 31):
        for i in range(500):
            sequences.append(_cycled_seq(i, length))
    return sequences


def _write_samplesheet(path, rows):
    lines = ["sample_id,fastq_1"]
    for sample_id, fastq_path in rows:
        lines.append(f"{sample_id},{fastq_path}")
    path.write_text("\n".join(lines) + "\n")


def test_samplesheet_processes_local_fastqs_without_download(tmp_path):
    fq1 = tmp_path / "s1.fastq"
    fq2 = tmp_path / "s2.fastq"
    _write_fastq(fq1, [_cycled_seq(i, 29) for i in range(200)])
    _write_fastq(fq2, [_cycled_seq(i, 30) for i in range(200)])

    sheet = tmp_path / "samplesheet.csv"
    _write_samplesheet(sheet, [("s1", fq1), ("s2", fq2)])

    out_dir = tmp_path / "out"
    result = run_cohort(sheet, out_dir, infer_reads=10_000)

    assert {e["sample_id"] for e in result.evidence} == {"s1", "s2"}
    for name, path in result.cohort_tsv_paths.items():
        assert path.exists()
    assert (out_dir / "s1.evidence.json").exists()
    assert (out_dir / "s2.evidence.json").exists()


def test_audit_mode_proposed_rules_match_production_applied_rules(tmp_path):
    sequences = _artifact_sequences()

    fq = tmp_path / "artifact.fastq"
    _write_fastq(fq, sequences)

    sheet = tmp_path / "samplesheet.csv"
    _write_samplesheet(sheet, [("artifact_sample", fq)])

    audit_out = tmp_path / "audit"
    audit_result = run_cohort(sheet, audit_out, audit_only=True, infer_reads=100_000)

    prod_out = tmp_path / "prod"
    prod_result = run_cohort(sheet, prod_out, audit_only=False, infer_reads=100_000)

    audit_evidence = audit_result.evidence[0]
    prod_evidence = prod_result.evidence[0]

    # Audit mode applied nothing.
    assert audit_evidence["trim_rules_applied"] == []
    proposed_by_length = {
        r["read_length_before"]: r["trim_bases"]
        for r in audit_evidence["trim_rules_proposed_not_applied"]
        if r.get("decision") == "safe_to_trim" or r.get("trim_bases")
    }

    applied_by_length = {
        r["read_length_before"]: r["trim_bases"]
        for r in prod_evidence["trim_rules_applied"]
    }

    assert applied_by_length == {30: 1}
    # What production actually applied, audit correctly identified as the
    # rule it would have applied (same length/trim_bases).
    assert proposed_by_length.get(30) == applied_by_length[30]

    # Audit mode really didn't modify reads: length-30 artifacts are still
    # length 30 in its output; production trimmed them to 29.
    audit_lengths = {
        len(line)
        for line in (audit_out / "artifact_sample.collapsed.fa").read_text().splitlines()
        if line and not line.startswith(">")
    }
    prod_lengths = {
        len(line)
        for line in (prod_out / "artifact_sample.collapsed.fa").read_text().splitlines()
        if line and not line.startswith(">")
    }
    assert 30 in audit_lengths
    assert 30 not in prod_lengths
