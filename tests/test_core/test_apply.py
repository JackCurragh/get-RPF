"""M3 acceptance tests: apply path + provenance.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M3: threading
per-length override_trims into RPFExtractor.extract_rpfs, and the M2->M3
boundary plan (core/apply.py) that drives it.
"""

from getRPF.core.apply import plan_from_file
from getRPF.core.processors.rpf_extractor import RPFExtractor

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def _write_fastq(path, sequences):
    lines = []
    for i, seq in enumerate(sequences):
        lines.append(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    path.write_text("".join(lines))


def _build_mixed_fixture(tmp_path):
    """20000 reads of length 30 with a 99% 3' 'A' artifact, plus small,
    clean control classes at 28/29/31 with no terminal signal."""
    sequences = []
    for i in range(20000):
        body = _cycled_seq(i, 29)
        tail = "A" if i < 19800 else BASES[i % 4]  # 99%
        sequences.append(body + tail)
    for length in (28, 29, 31):
        for i in range(500):
            sequences.append(_cycled_seq(i, length))

    input_file = tmp_path / "mixed.fastq"
    _write_fastq(input_file, sequences)
    return input_file


def _empty_architecture_extractor() -> RPFExtractor:
    extractor = RPFExtractor()
    extractor.architecture_db.architectures = []
    return extractor


def test_applied_rules_reproducible_from_provenance(tmp_path):
    input_file = _build_mixed_fixture(tmp_path)

    override_trims, applied_rules, proposed_rules, _decisions, _sketch = plan_from_file(
        input_file, format="fastq", infer_reads=100_000
    )

    assert override_trims == {30: {"trim_5p": 0, "trim_3p": 1}}
    assert len(applied_rules) == 1
    assert applied_rules[0]["read_length_before"] == 30
    assert applied_rules[0]["trim_bases"] == 1

    out1 = tmp_path / "out1.fastq"
    _empty_architecture_extractor().extract_rpfs(
        input_file,
        out1,
        format="fastq",
        collapsed_only=True,
        override_trims=override_trims,
    )

    # Reconstruct override_trims purely from the provenance records, as a
    # downstream consumer (e.g. --rules override mode) would.
    reconstructed = {
        r["read_length_before"]: {"trim_5p": 0, "trim_3p": r["trim_bases"]}
        for r in applied_rules
    }
    out2 = tmp_path / "out2.fastq"
    _empty_architecture_extractor().extract_rpfs(
        input_file,
        out2,
        format="fastq",
        collapsed_only=True,
        override_trims=reconstructed,
    )

    assert (tmp_path / "out1.collapsed.fa").read_text() == (
        tmp_path / "out2.collapsed.fa"
    ).read_text()


def test_before_after_artifact_removed_body_untouched(tmp_path):
    input_file = _build_mixed_fixture(tmp_path)

    override_trims, applied_rules, _proposed, _decisions, _sketch = plan_from_file(
        input_file, format="fastq", infer_reads=100_000
    )

    output_file = tmp_path / "out.fastq"
    _empty_architecture_extractor().extract_rpfs(
        input_file,
        output_file,
        format="fastq",
        collapsed_only=True,
        override_trims=override_trims,
    )

    collapsed_seqs = set()
    lines = (tmp_path / "out.collapsed.fa").read_text().splitlines()
    for i in range(0, len(lines), 2):
        collapsed_seqs.add(lines[i + 1])

    # The artifact class (length 30) must come out at length 29, with the
    # terminal 'A' removed, and body sequence untouched.
    expected_trimmed_bodies = {_cycled_seq(i, 29) for i in range(20000)}
    assert collapsed_seqs & expected_trimmed_bodies
    assert all(len(s) != 30 for s in collapsed_seqs)  # no un-trimmed 30-mers survive

    # The clean control classes (28/29/31) must pass through byte-identical
    # -- no trim was proposed for them.
    for length in (28, 29, 31):
        for i in range(500):
            assert _cycled_seq(i, length) in collapsed_seqs
