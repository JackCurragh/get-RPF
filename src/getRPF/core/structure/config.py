"""Every threshold used by structure inference, named and documented (spec §8).

Reports print the measured value next to each decision these control, so any
threshold can be checked against the data that crossed it.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class InferenceConfig:
    """Operational thresholds. None of them is a scientific constant."""

    sample_reads: int = 300_000
    """Reads in the bounded inference sample. Inference never needs more,
    because the architecture is shared by every read."""

    seed_k: int = 14
    """Length of the k-mer used to group reads that share a biological fragment."""

    top_seeds: int = 4000
    """How many of the most frequent k-mers may become group seeds."""

    max_seed_base_fraction: float = 0.6
    """A seed whose most common base exceeds this fraction is low complexity."""

    max_seed_homopolymer: int = 5
    """A seed with a longer homopolymer run is low complexity."""

    min_group_distinct: int = 5
    """Distinct sequences a group needs before it counts as informative."""

    min_others: int = 3
    """Other distinct sequences that must cover an offset before a base is scored."""

    min_groups: int = 100
    """Below this many informative groups the pileup is underpowered."""

    min_scored_per_position: int = 100
    """Scored bases a position needs before its agreement is reported."""

    window: int = 16
    """Positions examined inward from the read start and from the anchor."""

    agree_bio: float = 0.9
    """Within-fragment agreement at or above this is biology level."""

    agree_random: float = 0.4
    """Within-fragment agreement at or below this is random level (chance is 0.25)."""

    random_min_entropy: float = 1.7
    """Library-wide entropy (bits) a low-agreement position needs to count as
    random rather than skewed (non-templated-like)."""

    fixed_base_fraction: float = 0.9
    """A position whose most common base exceeds this fraction across the
    library is constant technical sequence, whatever its agreement."""

    composition_sample: int = 50_000
    """Reads used for library-wide base composition per position."""

    few_values_sample: int = 5_000
    """Reads used to test technical segments for a small set of values (barcodes)."""

    min_value_segment: int = 3
    """Shortest segment tested for few values (random 3-mers have 64 values)."""

    few_values_max: int = 12
    """A segment whose top values reaching ``few_values_coverage`` number at
    most this many is barcode-like."""

    few_values_coverage: float = 0.9
    """Fraction of reads the top values must cover for the few-values test."""

    few_values_growth: float = 1.5
    """A barcode segment is extended one base at a time while the number of
    values needed grows by no more than this factor (a random base multiplies
    it by about four)."""
