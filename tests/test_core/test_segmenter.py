"""Focused tests for the small HMM segmentation API."""

from getRPF.core.processors.segmenter import (
    ProbabilisticSegmenter,
    STATE_RPF,
    STATE_UMI,
)
from getRPF.core.processors.signals import SignalStats


def _stats(length=8):
    return SignalStats(
        entropy_5p=[1.8] * length,
        composition_5p=[{"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}] * length,
        dinucleotide_5p=[],
        entropy_3p=[1.8] * length,
        composition_3p=[{"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}] * length,
        dinucleotide_3p=[],
        sample_size=10,
    )


def test_viterbi_path_has_one_state_per_observation():
    segmenter = ProbabilisticSegmenter()

    path = segmenter._viterbi_path(segmenter._prepare_observations(_stats()))

    assert len(path) == len(_stats().entropy_5p)
    assert path[0] in {STATE_UMI, STATE_RPF}


def test_decode_with_posteriors_returns_aligned_viterbi_path():
    path, posteriors, segments = ProbabilisticSegmenter().decode_with_posteriors(
        _stats()
    )

    assert len(path) == len(posteriors) == len(_stats().entropy_5p)
    assert segments
    for position in posteriors:
        assert abs(sum(position) - 1.0) < 1e-9

