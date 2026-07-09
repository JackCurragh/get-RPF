"""M2 acceptance tests: boundary estimators, concordance, and refusal.

Covers docs/release_qc_and_terminal_trimming_plan.md Milestone M2 / section 6.
"""

from getRPF.core.processors.boundary import (
    _decide,
    concordance_and_confidence,
    estimate_boundaries,
)
from getRPF.core.processors.signals import SignalProcessor, SignalStats

BASES = "ACGT"


def _cycled_seq(i: int, length: int) -> str:
    return "".join(BASES[(i + pos) % 4] for pos in range(length))


def _per_length_stats(reads_by_length):
    proc = SignalProcessor()
    stats = {length: proc.process_reads(reads) for length, reads in reads_by_length.items()}
    support = {length: len(reads) for length, reads in reads_by_length.items()}
    return stats, support


def test_single_terminal_artifact_yields_exactly_one_trim():
    """A 3' single-base terminal artifact confined to one length class
    should trim exactly that class by exactly 1 base, and leave clean
    classes alone."""
    reads_by_length = {}

    # Clean classes: no terminal bias.
    for length in (28, 29, 31):
        reads_by_length[length] = [_cycled_seq(i, length) for i in range(1000)]

    # Artifact class: 99% of reads end in 'A'.
    artifact_reads = []
    for i in range(1000):
        body = _cycled_seq(i, 29)
        tail = "A" if i < 990 else BASES[i % 4]
        artifact_reads.append(body + tail)
    reads_by_length[30] = artifact_reads

    stats, support = _per_length_stats(reads_by_length)
    decisions = estimate_boundaries(
        stats, support, end="3p", s_min=1000, support_floor=200
    )

    assert decisions[30].decision == "safe_to_trim"
    assert decisions[30].trim_bases == 1
    assert decisions[30].n_indep >= 2

    for length in (28, 29, 31):
        assert decisions[length].decision != "safe_to_trim"


def test_internal_motif_proposes_nothing():
    """A repeated motif in the middle of the read (not touching either
    terminus) must not trigger a trim proposal on either end."""
    reads = []
    for i in range(1000):
        left = _cycled_seq(i, 10)
        middle = "AAAAAA"  # constant internal motif, positions 10-15
        right = _cycled_seq(i + 7, 10)
        reads.append(left + middle + right)

    stats, support = _per_length_stats({26: reads})

    for end in ("3p", "5p"):
        decisions = estimate_boundaries(
            stats, support, end=end, s_min=1000, support_floor=200
        )
        d = decisions[26]
        assert d.estimates == {}
        assert d.decision == "hold"


def test_divergent_estimators_do_not_trim_unit():
    """Direct test of the concordance/refusal core: when estimators
    disagree by a wide spread, confidence collapses and the decision must
    not be safe_to_trim, regardless of support."""
    estimates = {"b_kmer": 1, "b_entropy": 6}
    spread, n_indep, concordance, confidence = concordance_and_confidence(
        estimates, support=5000, s_min=1000, support_floor=200
    )
    assert spread == 5
    assert concordance == 0.0
    assert confidence == 0.0

    decision, *_ = _decide(
        estimates,
        support=5000,
        end="3p",
        named_element=False,
        tau_high_3p=0.7,
        tau_high_5p=0.85,
        s_min=1000,
        support_floor=200,
    )
    assert decision != "safe_to_trim"


def test_divergent_estimators_do_not_trim_end_to_end():
    """Same refusal behavior, exercised through estimate_boundaries with
    genuinely conflicting per-position signals (high dominant-base
    frequency out to position 4, but entropy only low at position 0)."""
    stats = SignalStats(
        entropy_5p=[],
        composition_5p=[],
        dinucleotide_5p=[],
        entropy_3p=[0.05] + [1.8] * 11,
        composition_3p=(
            [{"A": 0.9, "C": 0.033, "G": 0.033, "T": 0.034}] * 5
            + [{"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}] * 7
        ),
        dinucleotide_3p=[],
        sample_size=5000,
    )
    decisions = estimate_boundaries(
        {30: stats}, {30: 5000}, end="3p", s_min=1000, support_floor=200
    )
    d = decisions[30]
    assert d.estimates.get("b_kmer") == 5
    assert d.estimates.get("b_entropy") == 1
    assert d.spread == 4
    assert d.decision != "safe_to_trim"


def test_five_prime_never_trims_without_named_element():
    """Even with strong, concordant de novo 3'-style evidence, the 5' end
    must not reach safe_to_trim unless the length class has a named
    element (seqspec / high-posterior HMM segment)."""
    reads = []
    for i in range(20000):
        head = "A"  # strong 5' bias
        rest = _cycled_seq(i, 27)
        reads.append(head + rest)

    stats, support = _per_length_stats({28: reads})

    decisions_no_named = estimate_boundaries(
        stats, support, end="5p", s_min=1000, support_floor=200
    )
    assert decisions_no_named[28].decision != "safe_to_trim"

    decisions_named = estimate_boundaries(
        stats, support, end="5p", s_min=1000, support_floor=200,
        named_element_lengths={28},
    )
    assert decisions_named[28].decision == "safe_to_trim"
