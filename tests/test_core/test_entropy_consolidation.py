"""M0 acceptance tests: entropy consolidation onto SignalProcessor.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M0:
CleanlinessChecker's per-position checks must be coverage-normalized (reads
actually covering a position), not normalized by total sample reads, so
variable-length tails don't produce false low-complexity/entropy failures.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from getRPF.core.checkers import (
    InformationContentCheck,
    Status,
    categorize_failures,
    run_all_cleanliness_checks,
)
from getRPF.core.processors.check import CleanlinessChecker

BASES = "ACGT"


def _write_fastq(tmp_path, name, sequences):
    records = [
        SeqRecord(
            Seq(seq),
            id=f"read{i}",
            description="",
            letter_annotations={"phred_quality": [40] * len(seq)},
        )
        for i, seq in enumerate(sequences)
    ]
    path = tmp_path / name
    SeqIO.write(records, path, "fastq")
    return path


def _cycled_seq(i: int, length: int) -> str:
    """Deterministic sequence where every position cycles evenly through
    A/C/G/T as the read index i varies, so per-position composition among
    any contiguous block of reads is (close to) maximally diverse."""
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def test_low_coverage_tail_position_no_longer_fails(tmp_path):
    """A max-complexity tail position covered by only a few reads must not
    be flagged low-complexity just because it's rare in the whole sample."""
    sequences = []
    # 970 reads of length 28: the common body.
    for i in range(970):
        sequences.append(_cycled_seq(i, 28))
    # 30 reads of length 34: positions 28-31 are only covered by these reads,
    # but are still evenly cycled across A/C/G/T (maximal diversity).
    for i in range(970, 1000):
        sequences.append(_cycled_seq(i, 34))

    path = _write_fastq(tmp_path, "tail.fastq", sequences)
    checker = CleanlinessChecker(format="fastq", max_reads=None)
    results = checker.analyze_file(path)

    check_result = InformationContentCheck().check(results)
    low_positions = {pos for pos, _ in check_result.details["low_entropy_positions"]}

    assert not low_positions & {28, 29, 30, 31}
    assert check_result.status == Status.PASS


def test_mixed_length_clean_rpfs_pass(tmp_path):
    """Mixed 28-32nt reads with real per-position diversity should pass all
    cleanliness checks, including at the tail positions unique to 32-mers."""
    lengths = [28, 29, 30, 31, 32]
    sequences = [
        _cycled_seq(i, lengths[i % len(lengths)]) for i in range(1000)
    ]

    path = _write_fastq(tmp_path, "mixed.fastq", sequences)
    checker = CleanlinessChecker(format="fastq", max_reads=None)
    results = checker.analyze_file(path)

    check_results = run_all_cleanliness_checks(results)
    failing = {name: r.message for name, r in check_results.items() if r.status == Status.FAIL}

    assert failing == {}
    assert categorize_failures(check_results)["is_clean"] is True


def test_synthetic_three_prime_terminal_bias_flagged(tmp_path):
    """Reads with a dominant 3' terminal base should be flagged as
    three_prime_terminal_bias, distinct from generic end_bias."""
    sequences = []
    for i in range(1000):
        body = _cycled_seq(i, 28)
        # 99% of reads end in 'A'; the rest cycle through other bases.
        tail = "A" if i < 990 else BASES[i % 4]
        sequences.append(body + tail)

    path = _write_fastq(tmp_path, "terminal_bias.fastq", sequences)
    checker = CleanlinessChecker(format="fastq", max_reads=None)
    results = checker.analyze_file(path)

    check_results = run_all_cleanliness_checks(results)
    assert check_results["end_bias"].status == Status.FAIL
    assert check_results["end_bias"].details["worst_position"].startswith("3'")

    categories = categorize_failures(check_results)
    assert "three_prime_terminal_bias" in categories["failure_categories"]
    assert "five_prime_terminal_bias" not in categories["failure_categories"]
    assert "end_bias" not in categories["failure_categories"]
