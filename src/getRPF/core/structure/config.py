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

    # --- Steps 0-1: observation and profile (spec §5.0-5.1) ---

    raw_modal_fraction: float = 0.9
    """At least this fraction of reads at one length: raw fixed-length cycles."""

    trimmed_modal_fraction: float = 0.5
    """At most this fraction at the modal length: already trimmed."""

    # --- Step 2: anchor (spec §5.2) ---

    anchor_seed: int = 10
    """Exact seed length used to find adapter candidates. Two seeds are tried
    (adapter positions 1-10 and 11-20), so one error in either is tolerated."""

    min_partial_overlap: int = 10
    """Shortest adapter overlap accepted where the adapter runs off the read end."""

    max_mismatch_per_10nt: float = 1.0
    """Mismatches allowed per 10 nt of adapter overlap."""

    min_insert: int = 15
    """An adapter starting before this position marks a dimer, not an anchor."""

    min_anchor_support: float = 0.3
    """Fraction of reads an adapter must be located in to become the anchor."""

    max_index_length: int = 12
    """Longest index expected between the adapter and P7 (chain check)."""

    tail_min_run: int = 6
    """Homopolymer run directly before the adapter that counts as a tail."""

    tail_min_fraction: float = 0.5
    """Fraction of anchored reads that must carry such a run for a tail
    architecture to be declared."""

    typical_footprint: int = 30
    """Typical monosome footprint, used only to size search windows (spec §2
    rule 3). It never places a boundary."""

    max_window: int = 40
    """Upper bound on the search window derived from the anchor arithmetic."""

    denovo_sample: int = 50_000
    """Reads used for de novo adapter discovery."""

    denovo_k: int = 12
    """k-mer length for de novo adapter discovery."""

    denovo_min_support: float = 0.1
    """Fraction of reads a de novo candidate k-mer must occur in."""

    denovo_min_position_sd: float = 0.75
    """A de novo candidate's position must vary this much (nt, SD) across
    reads. An adapter's position varies with insert length (SD ~1.4 nt for
    28-32 nt footprints); a fixed 5' element's barely varies at all."""

    denovo_extend_dominance: float = 0.8
    """A de novo element is extended while the next base agrees in this
    fraction of occurrences."""

    denovo_max_length: int = 34
    """Longest de novo adapter reported."""
