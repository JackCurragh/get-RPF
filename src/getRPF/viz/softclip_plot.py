"""Soft-clipping visualization for alignment/STAR mode.

Two display modes:
1) Summary (JSON only): per-read-length mean soft clips (5'/3') and recommended trims.
2) Heatmap (JSON + BAM): for each read length, distribution of 5' and 3' clip sizes
   (fraction of reads at each clip value), overlaid with recommended trims.

Use with CLI: getRPF plot-softclips --align-json <json> [--bam <bam>] -o <png>
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import pysam  # optional for heatmap mode

    PYSAM_AVAILABLE = True
except Exception:
    PYSAM_AVAILABLE = False


def _load_align_json(path: Path) -> Dict[str, Any]:
    """Load either align-detect JSON OR extraction_report JSON.

    Normalizes keys to: stats (alignment_statistics), trims (recommended_*), perlen (per-length dict).
    """
    data = json.loads(Path(path).read_text())

    # Case 1: align-detect JSON
    if "per_length_analysis" in data or "trim_recommendations" in data:
        stats = data.get("alignment_statistics", {})
        trims = data.get("trim_recommendations", {})
        perlen = data.get("per_length_analysis", {})
        return {"stats": stats, "trims": trims, "perlen": perlen}

    # Case 2: extraction_report JSON from alignment-based extractor
    if "trim_boundaries" in data or "extraction_summary" in data:
        stats = data.get("alignment_statistics", {})
        tb = data.get("trim_boundaries", {})
        # Map consensus trims to the same fields used by the plotter
        conf_map = {"high": 0.9, "medium": 0.7, "low": 0.5}
        trims = {
            "recommended_5prime_trim": tb.get("consensus_5p", 0),
            "recommended_3prime_trim": tb.get("consensus_3p", 0),
            "consensus_level": conf_map.get(tb.get("confidence")),
        }
        # Per-length shape differs: use what we have; means may be missing
        perlen_raw = tb.get("per_length", {}) or {}
        perlen = {}
        for k, row in perlen_raw.items():
            # Ensure string keys and provide placeholders for means
            perlen[str(k)] = {
                "mean_5prime_clips": row.get("trim_5p", 0),  # proxy
                "mean_3prime_clips": row.get("trim_3p", 0),  # proxy
                "5prime_trim": row.get("trim_5p", 0),
                "3prime_trim": row.get("trim_3p", 0),
                "n_reads": row.get("n_reads", 0),
            }
        return {"stats": stats, "trims": trims, "perlen": perlen}

    # Fallback empty
    return {"stats": {}, "trims": {}, "perlen": {}}


def _compute_clip_distributions_from_bam(bam_path: Path):
    if not PYSAM_AVAILABLE:
        raise RuntimeError(
            "pysam is required for heatmap mode; install pysam or disable --heatmap."
        )

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        # Determine bounds for read length and max clip
        lengths = []
        max_clip5 = 0
        max_clip3 = 0
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or not read.cigartuples:
                continue
            ln = read.query_length
            lengths.append(ln)
            cig = read.cigartuples
            c5 = cig[0][1] if cig[0][0] == 4 else 0
            c3 = cig[-1][1] if cig[-1][0] == 4 else 0
            max_clip5 = max(max_clip5, c5)
            max_clip3 = max(max_clip3, c3)

        if not lengths:
            return None

        min_len, max_len = min(lengths), max(lengths)
        # Re-scan to fill histograms
        heat5 = np.zeros((max_len - min_len + 1, max_clip5 + 1), dtype=np.int64)
        heat3 = np.zeros((max_len - min_len + 1, max_clip3 + 1), dtype=np.int64)
        counts = np.zeros((max_len - min_len + 1,), dtype=np.int64)

        bam.reset()
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or not read.cigartuples:
                continue
            ln = read.query_length
            idx = ln - min_len
            cig = read.cigartuples
            c5 = cig[0][1] if cig[0][0] == 4 else 0
            c3 = cig[-1][1] if cig[-1][0] == 4 else 0
            heat5[idx, c5] += 1
            heat3[idx, c3] += 1
            counts[idx] += 1

        # Convert to fractions per row
        with np.errstate(divide="ignore", invalid="ignore"):
            heat5 = heat5 / counts[:, None]
            heat3 = heat3 / counts[:, None]
            heat5 = np.nan_to_num(heat5)
            heat3 = np.nan_to_num(heat3)

        return {
            "min_len": int(min_len),
            "max_len": int(max_len),
            "heat5": heat5,
            "heat3": heat3,
        }


def plot_softclips(
    align_json: Path,
    output_png: Path,
    bam_path: Optional[Path] = None,
    heatmap: bool = True,
    title: Optional[str] = None,
):
    jd = _load_align_json(Path(align_json))
    trims, perlen = jd["trims"], jd["perlen"]

    # Prepare per-length arrays from summary
    lengths = []
    mean5 = []
    mean3 = []
    rec5 = []
    rec3 = []
    for Ls, row in perlen.items():
        L = int(Ls)
        lengths.append(L)
        mean5.append(row.get("mean_5prime_clips", 0))
        mean3.append(row.get("mean_3prime_clips", 0))
        rec5.append(row.get("5prime_trim", 0))
        rec3.append(row.get("3prime_trim", 0))

    # Sort by length for nice plotting
    idx = np.argsort(lengths) if lengths else []
    lengths = np.array(lengths)[idx] if len(idx) else np.array([])
    mean5 = np.array(mean5)[idx] if len(idx) else np.array([])
    mean3 = np.array(mean3)[idx] if len(idx) else np.array([])
    rec5 = np.array(rec5)[idx] if len(idx) else np.array([])
    rec3 = np.array(rec3)[idx] if len(idx) else np.array([])

    # If no per-length JSON but we do have BAM, render heatmap-only layout
    has_heatmap = heatmap and (bam_path is not None)
    if lengths.size == 0 and has_heatmap:
        H = _compute_clip_distributions_from_bam(Path(bam_path))
        fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
        if not H:
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                "No aligned reads in BAM for heatmap",
                ha="center",
                va="center",
            )
        else:
            min_len, max_len = H["min_len"], H["max_len"]
            heat5 = H["heat5"]
            im = ax.imshow(
                heat5,
                aspect="auto",
                origin="lower",
                extent=[0, heat5.shape[1], min_len, max_len],
                cmap="viridis",
            )
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label("Fraction of reads")
            ax.set_xlabel("5' soft clip (bases)")
            ax.set_ylabel("Read length")
            # Overlay global 5' trim line if within range
            g5 = trims.get("recommended_5prime_trim", 0)
            if g5 <= max(1, heat5.shape[1] - 1):
                ax.axvline(g5, color="w", ls="--", lw=1.2, alpha=0.8)
        if title:
            ax.set_title(title, pad=10)
        output_png = Path(output_png)
        output_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_png, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return output_png

    # Layout when per-length arrays are present
    ncols = 3 if has_heatmap else 2
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 4), constrained_layout=True)
    if ncols == 2:
        ax0, ax1 = axes
    else:
        ax0, ax1, ax2 = axes

    # Panel 1: per-length means
    if lengths.size:
        ax0.plot(lengths, mean5, label="Mean 5' clips", color="#1f77b4", lw=1.8)
        ax0.plot(lengths, mean3, label="Mean 3' clips", color="#ff7f0e", lw=1.8)
        ax0.set_xlabel("Read length")
        ax0.set_ylabel("Mean soft clips (bases)")
        ax0.legend(frameon=False)
    ax0.set_title("Per-length mean soft clips")

    # Panel 2: recommended trims per length
    if lengths.size:
        ax1.plot(lengths, rec5, label="5' trim", color="#1f77b4", lw=1.8)
        ax1.plot(lengths, rec3, label="3' trim", color="#ff7f0e", lw=1.8)
        # Global recommended trims (dashed)
        g5 = trims.get("recommended_5prime_trim", 0)
        g3 = trims.get("recommended_3prime_trim", 0)
        ax1.axhline(g5, color="#1f77b4", ls="--", alpha=0.6)
        ax1.axhline(g3, color="#ff7f0e", ls="--", alpha=0.6)
        ax1.text(
            0.02,
            0.95,
            f"global 5'={g5}, 3'={g3}\nconsensus={trims.get('consensus_level',0):.2f}",
            transform=ax1.transAxes,
            va="top",
            ha="left",
            fontsize=9,
        )
        ax1.set_xlabel("Read length")
        ax1.set_ylabel("Recommended trim (bases)")
        ax1.legend(frameon=False)
    ax1.set_title("Per-length recommended trims")

    # Panel 3: heatmaps from BAM (optional)
    if has_heatmap:
        try:
            H = _compute_clip_distributions_from_bam(Path(bam_path))
        except Exception as e:
            ax2.axis("off")
            ax2.text(0.5, 0.5, f"Heatmap disabled: {e}", ha="center", va="center")
        else:
            if not H:
                ax2.axis("off")
                ax2.text(
                    0.5, 0.5, "No aligned reads for heatmap", ha="center", va="center"
                )
            else:
                min_len, max_len = H["min_len"], H["max_len"]
                heat5 = H["heat5"]
                # Build a 2-row heatmap (5' top, 3' bottom) stacked vertically in same axis
                # Simpler: show 5' only; 3' only; or combine—as a quick MVP show 5'.
                im = ax2.imshow(
                    heat5,
                    aspect="auto",
                    origin="lower",
                    extent=[0, heat5.shape[1], min_len, max_len],
                    cmap="viridis",
                )
                cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
                cb.set_label("Fraction of reads")
                ax2.set_xlabel("5' soft clip (bases)")
                ax2.set_ylabel("Read length")
                ax2.set_title("5' clip distribution heatmap")
                # Overlay global 5' trim line if within range
                g5 = trims.get("recommended_5prime_trim", 0)
                if g5 <= heat5.shape[1]:
                    ax2.axvline(g5, color="w", ls="--", lw=1.2, alpha=0.8)

    if title:
        fig.suptitle(title, y=1.02)

    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_png
