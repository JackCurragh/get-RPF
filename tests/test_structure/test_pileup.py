"""P3 gate: the pileup must recover known structure, and no more than that.

Anchors come from the simulation truth, so these tests exercise the pileup
alone (spec §11). A test fails if the answer is wrong *or* more certain than
the truth allows.
"""

from collections import Counter

from getRPF.core.structure import Status, infer_junctions
from getRPF.core.structure.config import InferenceConfig
from getRPF.core.structure.pileup import _score, _Tally

from .sim import barcode, fixed, nta, random_bases, simulate


def infer(library, trimmed=False):
    return infer_junctions(library.reads, None if trimmed else library.anchors)


def contradicted(answer):
    return {alt.value for alt in answer.alternatives if alt.verdict == "contradicted"}


def test_case1_insert_then_adapter_has_no_technical_bases():
    result = infer(simulate())
    assert result.pileup.groups >= 100
    for answer in (result.q2, result.q3):
        assert answer.status is Status.RESOLVED
        assert answer.value.technical_length == 0
        assert answer.value.nta_length == 0
        assert 1 in contradicted(answer)


def test_case2_five_prime_umi():
    result = infer(simulate(five_prime=[random_bases(5)]))
    assert result.q2.status is Status.RESOLVED
    assert result.q2.value.technical_length == 5
    assert result.q2.value.blocks == (("random", 5),)
    assert {4, 6} <= contradicted(result.q2)
    assert result.q3.value.technical_length == 0


def test_case3_three_prime_random_block():
    result = infer(simulate(three_prime=[random_bases(4)]))
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.technical_length == 4
    assert result.q3.value.blocks == (("random", 4),)
    assert result.q2.value.technical_length == 0


def test_case4_random_blocks_in_both_frames():
    result = infer(
        simulate(five_prime=[random_bases(3)], three_prime=[random_bases(4)])
    )
    assert (result.q2.status, result.q3.status) == (Status.RESOLVED, Status.RESOLVED)
    assert result.q2.value.technical_length == 3
    assert result.q3.value.technical_length == 4
    assert result.first_pass.consensus_trim == (0, 0)
    assert result.pileup.consensus_trim == (3, 4)


def test_case6_five_prime_nta_is_an_interval_not_a_trim():
    skewed_g = {"G": 0.85, "A": 0.05, "C": 0.05, "T": 0.05}
    result = infer(simulate(five_prime=[nta({0: 0.3, 1: 0.5, 2: 0.2}, skewed_g)]))
    call = result.q2.value
    assert result.q2.status is Status.INTERVAL
    assert call.technical_length == 0
    assert 1 <= call.nta_length <= 2
    assert max(call.nta_bases, key=call.nta_bases.get) == "G"


def test_case7_three_prime_nta():
    skewed_a = {"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1}
    result = infer(simulate(three_prime=[nta({1: 1.0}, skewed_a)]))
    call = result.q3.value
    assert result.q3.status is Status.INTERVAL
    assert call.technical_length == 0
    assert call.nta_length == 1


def test_case8_random_linker_then_single_barcode():
    result = infer(simulate(three_prime=[random_bases(5), fixed("TGCAT")]))
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.technical_length == 10
    assert result.q3.value.blocks == (("fixed", 5), ("random", 5))


def test_case9_random_linker_then_multiplexed_barcodes():
    four = barcode(["ACGTA", "CGTAC", "GTACG", "TACGT"])
    result = infer(simulate(three_prime=[random_bases(5), four]))
    assert result.q3.status is Status.RESOLVED
    assert result.q3.value.technical_length == 10
    assert result.q3.value.blocks == (("barcode", 5), ("random", 5))


def test_case10_trimmed_reads_use_the_read_end_as_anchor():
    result = infer(
        simulate(five_prime=[random_bases(4)], read_length=None), trimmed=True
    )
    assert result.q2.value.technical_length == 4
    assert result.q3.value.technical_length == 0


def test_case13_low_duplication_is_underpowered_not_absent():
    library = simulate(
        five_prime=[random_bases(5)], n_reads=20_000, n_transcripts=5_000, zipf_s=0.0
    )
    result = infer(library)
    for answer in (result.q2, result.q3):
        assert answer.status is Status.UNDERPOWERED
        assert answer.value is None
        assert "absent" not in answer.explanation


def test_five_prime_profile_uses_reads_without_a_located_three_prime_anchor():
    library = simulate(five_prime=[random_bases(5)])
    config = InferenceConfig(
        min_groups=3, min_group_distinct=2, min_scored_per_position=3
    )
    result = infer_junctions(library.reads, [None] * len(library.reads), config)
    assert result.first_pass.reads_used == len(library.reads)
    assert result.first_pass.anchor_reads_used == 0
    assert result.q2.evidence[0].value > 0
    assert result.q3.status is Status.NOT_OBSERVABLE


def test_every_answer_carries_provenance():
    result = infer(simulate(five_prime=[random_bases(5)]))
    for answer in (result.q2, result.q3):
        assert answer.evidence and answer.alternatives
        assert all(alt.reason for alt in answer.alternatives)
        assert "groups" in answer.explanation
    assert (
        "5 technical bases from the read start (5 nt random)" in result.q2.explanation
    )


def test_tied_leave_one_out_consensus_gets_fractional_credit():
    tally = _Tally()
    _score(Counter({"A": 2, "C": 2}), "A", False, tally, min_others=3)
    assert tally.scored == 1
    assert tally.agree == 0.5


def test_four_way_tied_consensus_does_not_pass_random_cutoff():
    tally = _Tally()
    _score(Counter({"A": 1, "C": 1, "G": 1, "T": 1}), "A", False, tally, min_others=3)
    assert tally.scored == 1
    assert tally.agree == 0.25
