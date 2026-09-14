"""Alignment check logic on hand-made SAM records (no STAR needed)."""

from pathlib import Path

from getRPF.core.structure.align import (
    AlignmentProfile,
    _star_input_reads,
    check_junction,
    parse_sam,
    run_star,
)
from getRPF.core.structure.config import InferenceConfig
from getRPF.core.structure.model import Answer, JunctionCall, Status


def sam(flag, cigar, seq, md):
    return f"r\t{flag}\tchrI\t100\t255\t{cigar}\t*\t0\t0\t{seq}\t*\tNH:i:1\tMD:Z:{md}"


def test_parse_sam_orients_clips_and_bases_by_strand():
    lines = [
        "@HD\tVN:1.6",
        sam(0, "1S29M", "T" + "C" * 29, "29"),  # forward: 1 nt 5' clip, base T
        sam(16, "29M1S", "G" * 29 + "A", "10A18"),  # reverse: 5' clip is on the right
        sam(4, "*", "ACGT", ""),  # unmapped: ignored
        sam(256, "30M", "C" * 30, "30"),  # secondary: ignored
    ]
    profile = parse_sam(lines)
    assert profile.aligned == 2
    assert profile.clip5 == {1: 2}
    assert profile.clip3 == {0: 2}
    assert profile.clipped5_bases == {"T": 2}  # revcomp of the reverse read's last A
    assert profile.mismatches == 1 and profile.aligned_bases == 58


def answer(nta_length, nta_bases, status=Status.INTERVAL):
    call = JunctionCall("read_start", 0, (), nta_length, (0.4,) * nta_length, nta_bases)
    return Answer("Q2", call, status, (), (), "Pileup explanation.")


def profile(one, two, bases, aligned=1000, mismatches=5):
    p = AlignmentProfile(
        aligned=aligned, aligned_bases=30 * aligned, mismatches=mismatches
    )
    p.clip5.update({0: aligned - one - two, 1: one, 3: two})
    p.clipped5_bases.update(bases)
    return p


def test_junction_base_confirmed_when_clips_match_the_pileup():
    check = check_junction(
        answer(1, {"T": 0.6, "A": 0.4}),
        "5",
        profile(400, 0, {"T": 300, "A": 100}),
        InferenceConfig(),
    )
    assert check.verdict == "confirmed"
    assert check.answer.status is Status.INTERVAL
    assert "Alignment: alignment confirms" in check.answer.explanation


def test_junction_base_disagreement_is_flagged_not_blocking():
    check = check_junction(
        answer(1, {"T": 1.0}), "5", profile(5, 0, {"T": 5}), InferenceConfig()
    )
    assert check.verdict == "conflicting"
    assert check.answer.status is Status.INTERVAL
    assert check.flag is not None and check.flag.startswith("Q2:")


def test_missed_technical_block_is_a_conflict():
    check = check_junction(
        answer(0, {}, Status.RESOLVED),
        "5",
        profile(0, 500, {"A": 500}),
        InferenceConfig(),
    )
    assert check.verdict == "conflicting"
    assert check.answer.status is Status.CONFLICTING
    assert check.flag is None
    assert check.longer_clip_mode == 3


def test_run_star_with_nothing_to_align_needs_no_star():
    assert run_star([], Path("/nonexistent")).aligned == 0


def test_star_input_count_is_read_from_the_final_log(tmp_path):
    log = tmp_path / "Log.final.out"
    log.write_text("                          Number of input reads |\t0\n")
    assert _star_input_reads(log) == 0
    assert _star_input_reads(tmp_path / "missing.out") is None


def test_no_alignments_is_underpowered_not_a_contradiction():
    empty = AlignmentProfile(sequences=10_000)
    check = check_junction(answer(1, {"T": 1.0}), "5", empty, InferenceConfig())
    assert check.verdict == "underpowered"
    assert check.answer.status is Status.INTERVAL


def test_clean_junction_stays_resolved():
    check = check_junction(
        answer(0, {}, Status.RESOLVED),
        "5",
        profile(10, 5, {"A": 15}),
        InferenceConfig(),
    )
    assert check.verdict == "clean"
    assert check.answer.status is Status.RESOLVED
