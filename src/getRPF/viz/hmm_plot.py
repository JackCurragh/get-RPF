"""HMM entropy plot utility for getRPF.

Generates a PNG figure showing per-position Shannon entropy (5') with
color-coded HMM segments (UMI/RPF/ADAPTER/BARCODE when present).

Usage (via CLI wrapper):
  python3 -m getRPF.cli plot-hmm \
    input.fastq -f fastq -o hmm_plot.png --max-reads 20000

The function uses existing metric/segmenter utilities in getRPF:
  - SignalProcessor to compute per-position entropy and composition
  - ProbabilisticSegmenter (Viterbi) to infer segment boundaries
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib

# Use non-interactive backend for headless environments
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..core.processors.collapsed import TwoStageCollapser
from ..core.processors.segmenter import (
    decode_segments_with_posteriors,
    segment_reads,
)
from ..core.processors.signals import process_reads


def _sample_reads(
    input_file: Path,
    format: str,
    max_reads: int = 20000,
) -> List[str]:
    """Return up to `max_reads` sequences as uppercase strings.

    - For fastq/fasta: uses TwoStageCollapser.collapse_raw to stream and count
      unique reads, then expands proportionally up to `max_reads`.
    - For collapsed: respects header counts when sampling.
    """
    collapser = TwoStageCollapser()
    # Limit read scan to avoid pulling entire files for big runs
    counts = collapser.collapse_raw(
        Path(input_file), format=format, max_reads=max_reads
    )

    if not counts:
        return []

    total = sum(counts.values())
    target = min(max_reads, total)

    # Build a weighted reservoir sample without expanding huge counts in memory
    # Convert counts to a list of (seq, weight)
    items: List[Tuple[str, int]] = list(counts.items())

    # If total <= target, expand exactly
    if total <= target:
        reads: List[str] = []
        for seq, ct in items:
            reads.extend([seq] * ct)
        return [s[:512] for s in reads]  # cap length defensively

    # Otherwise, sample with probability proportional to count
    reads = []
    # Precompute cumulative weights for efficient sampling
    cumulative = []
    cum = 0
    for _, ct in items:
        cum += ct
        cumulative.append(cum)

    for _ in range(target):
        r = random.randint(1, total)
        # Binary search
        lo, hi = 0, len(cumulative) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cumulative[mid] < r:
                lo = mid + 1
            else:
                hi = mid
        reads.append(items[lo][0])

    return [s[:512] for s in reads]


def plot_hmm_entropy(
    input_file: Path,
    format: str,
    output_png: Path,
    max_reads: int = 20000,
    title: Optional[str] = None,
    show_segments: bool = False,
    show_freq: bool = True,
    show_posteriors: bool = False,
) -> Path:
    """Compute stats, run HMM segmentation, and save an entropy+segments PNG.

    Args:
        input_file: FASTQ/FASTA/collapsed path
        format: one of {'fastq','fasta','collapsed'}
        output_png: destination PNG path
        max_reads: cap on reads to sample for stats
        title: optional title for the figure

    Returns:
        Path to the written PNG.
    """
    reads = _sample_reads(Path(input_file), format=format, max_reads=max_reads)
    if not reads:
        raise RuntimeError("No reads available to plot HMM entropy.")

    # Compute per-position entropy/composition
    stats = process_reads(reads, compute_dinucleotide=False)

    # Segment via HMM (can be hidden with show_segments=False)
    if show_posteriors:
        _, posteriors, segments = decode_segments_with_posteriors(stats)
    else:
        posteriors, segments = [], segment_reads(stats)

    # Prepare plot
    x = list(range(len(stats.entropy_5p)))
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(x, stats.entropy_5p, color="#333333", lw=1.8, label="Entropy (5')")
    ax.set_xlabel("Read position (5' aligned)")
    ax.set_ylabel("Entropy (bits)")
    ax.set_ylim(0, max(2.1, max(stats.entropy_5p) if stats.entropy_5p else 2.0))
    if title:
        ax.set_title(title)

    lines = [plt.Line2D([], [], color="#333", lw=1.8, label="Entropy")]

    # Optional overlay: nucleotide frequencies on twin y-axis
    if show_freq and stats.composition_5p:
        ax2 = ax.twinx()
        ax2.set_ylabel("Nucleotide frequency")
        ax2.set_ylim(0.0, 1.0)
        # Build per-base series
        bases = ["A", "C", "G", "T"]
        colors_nt = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}
        nt_lines = []
        for b in bases:
            yb = [pos.get(b, 0.0) for pos in stats.composition_5p]
            (ln,) = ax2.plot(x, yb, lw=1.2, alpha=0.7, color=colors_nt[b], label=b)
            nt_lines.append(ln)
        lines.extend(nt_lines)

    # Optional overlay: HMM segments as translucent spans
    if show_segments and segments:
        colors = {
            "umi": "#f5a623",
            "barcode": "#bd10e0",
            "rpf": "#4a90e2",
            "adapter": "#d0021b",
        }
        for s in segments:
            c = colors.get(s.segment_type.lower(), "#999999")
            ax.axvspan(s.start_pos, s.end_pos, color=c, alpha=0.18, lw=0)

    # Optional: per-state posterior ribbons (stacked area)
    if show_posteriors and posteriors:
        import numpy as np

        P = np.array(posteriors)
        # We’ll plot UMI, RPF, ADAPTER; ignore START/END/BARCODE for clarity
        state_ix = {"UMI": 1, "RPF": 3, "ADAPTER": 4}
        cols = ["#f5a623", "#4a90e2", "#d0021b"]
        ys = [P[:, state_ix["UMI"]], P[:, state_ix["RPF"]], P[:, state_ix["ADAPTER"]]]
        ax.fill_between(x, 0, ys[0], color=cols[0], alpha=0.10, linewidth=0)
        ax.fill_between(x, ys[0], ys[0] + ys[1], color=cols[1], alpha=0.10, linewidth=0)
        ax.fill_between(
            x,
            ys[0] + ys[1],
            ys[0] + ys[1] + ys[2],
            color=cols[2],
            alpha=0.10,
            linewidth=0,
        )

    # Unified legend (entropy + bases; no segment legend by default)
    labs = [str(ln.get_label()) for ln in lines]
    ax.legend(lines, labs, loc="upper right", frameon=False, ncol=3)

    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png, dpi=200)
    plt.close(fig)
    return output_png
