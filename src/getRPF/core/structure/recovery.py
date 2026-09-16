"""Audit-only evidence for fixed-length read recovery (Phase 2 R0--R2).

This module deliberately never returns an approved production architecture.
It produces comparable candidates so a family decision can be reviewed and
then passed to the same :class:`Transform` used in production.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, Iterable, Sequence

from .anchors import CATALOGUE, probe_short_end
from .config import InferenceConfig


def _entropy(values: Iterable[str]) -> float:
    counts = Counter(values)
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts.values())


def audit_fixed_length_reads(
    reads: Sequence[str],
    qualities: Sequence[str] = (),
    config: InferenceConfig | None = None,
) -> Dict[str, Any]:
    """Return structural evidence and the predeclared trim candidate grid.

    Candidate trims are observations only. No candidate is selected from a
    modal insert length, and the short adapter probe remains separate from
    ordinary anchor inference.
    """
    config = config or InferenceConfig()
    if not reads:
        return {"reads": 0, "candidate_grid": [], "short_adapter_candidates": []}
    max_len = max(map(len, reads))
    lengths = sorted(set(map(len, reads)))
    candidates = []
    for left in range(13):
        for right in range(13):
            retained = [read[left : len(read) - right if right else None] for read in reads if len(read) > left + right]
            candidates.append(
                {
                    "left_trim": left,
                    "right_trim": right,
                    "n": len(retained),
                    "lengths": dict(sorted(Counter(map(len, retained)).items())),
                    "terminal_5p_kmer": Counter(read[:5] for read in retained).most_common(3),
                    "terminal_3p_kmer": Counter(read[-5:] for read in retained).most_common(3),
                }
            )
    quality_mean = None
    if qualities:
        scores = [ord(base) - 33 for quality in qualities for base in quality]
        quality_mean = round(sum(scores) / len(scores), 3) if scores else None
    return {
        "reads": len(reads),
        "read_lengths": lengths,
        "max_read_length": max_len,
        "base_composition_by_position": [dict(Counter(read[pos] for read in reads)) for pos in range(max_len)],
        "entropy_by_position": [_entropy(read[pos] for read in reads if len(read) > pos) for pos in range(max_len)],
        "quality_mean_phred": quality_mean,
        "short_adapter_candidates": [
            {
                "adapter": item.adapter.name,
                "sequence": item.adapter.sequence,
                "overlap": item.overlap,
                "support": round(item.support, 4),
                "mismatches": item.mismatches,
            }
            for item in probe_short_end(reads, config, CATALOGUE)
        ],
        "candidate_grid": candidates,
        "selection": "unresolved_audit_only",
    }
