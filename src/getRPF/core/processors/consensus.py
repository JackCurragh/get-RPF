"""Consensus processor for determining optimal trim parameters.

This module implements the logic to synthesize results from:
1. Architecture Detection (SeqSpec/Pattern Matching)
2. Alignment Verification (STAR Soft-clipping)

To produce a single, trustworthy set of trimming parameters.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)


@dataclass
class TrimConsensus:
    """Consensus trim parameters with confidence/reasoning."""

    trim_5p: int
    trim_3p: int
    adapter_sequence: Optional[str]
    confidence: float
    method: str  # "agreement", "architecture_dominant", "alignment_dominant"
    details: Dict[str, Any]


class TrimDecider:
    """Compatibility wrapper for :func:`decide_trim_consensus`."""

    def decide(
        self, architecture_result: Dict[str, Any], alignment_result: Dict[str, Any]
    ) -> TrimConsensus:
        return decide_trim_consensus(architecture_result, alignment_result)


def decide_trim_consensus(
    architecture_result: Dict[str, Any],
    alignment_result: Dict[str, Any],
) -> TrimConsensus:
    """Derive trim parameters from architecture and alignment evidence.

    This is intentionally a function: the decision has no mutable state and
    depends only on its two evidence objects. ``TrimDecider`` remains above as
    a compatibility wrapper for callers using the original API.
    """
    architecture_recommendations = architecture_result.get("trim_recommendations", {})
    architecture_5p = architecture_recommendations.get("recommended_5prime_trim", 0)
    architecture_3p = architecture_recommendations.get("recommended_3prime_trim", 0)
    architecture_adapter = architecture_recommendations.get("three_prime_adapter")
    architecture_match = architecture_result.get("architecture_match")

    alignment_recommendations = alignment_result.get("trim_recommendations", {})
    alignment_5p = alignment_recommendations.get("recommended_5prime_trim", 0)
    alignment_3p = alignment_recommendations.get("recommended_3prime_trim", 0)
    alignment_consensus = alignment_recommendations.get("consensus_level", 0.0)
    global_alignment_pattern = alignment_recommendations.get(
        "global_pattern_detected", False
    )

    if architecture_5p == alignment_5p:
        final_5p = architecture_5p
        method_5p = "agreement"
    elif architecture_match and architecture_5p > 0:
        final_5p = architecture_5p
        method_5p = "architecture_dominant"
    elif global_alignment_pattern and alignment_consensus > 0.8:
        final_5p = alignment_5p
        method_5p = "alignment_dominant"
    else:
        final_5p = architecture_5p if architecture_match else alignment_5p
        method_5p = "fallback"

    details = {
        "sources": {
            "architecture": {
                "5p": architecture_5p,
                "3p": architecture_3p,
                "match": architecture_match,
            },
            "alignment": {
                "5p": alignment_5p,
                "3p": alignment_3p,
                "consensus": alignment_consensus,
            },
        },
        "decisions": {"5p_logic": method_5p},
    }

    return TrimConsensus(
        trim_5p=final_5p,
        trim_3p=0,
        adapter_sequence=architecture_adapter,
        confidence=0.9 if method_5p == "agreement" else 0.7,
        method="combined_consensus",
        details=details,
    )
