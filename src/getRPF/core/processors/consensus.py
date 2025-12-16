"""Consensus processor for determining optimal trim parameters.

This module implements the logic to synthesize results from:
1. Architecture Detection (SeqSpec/Pattern Matching)
2. Alignment Verification (STAR Soft-clipping)

To produce a single, trustworthy set of trimming parameters.
"""

import logging
from typing import Dict, Any, Optional
from dataclasses import dataclass

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
    """Decides on trim parameters by comparing architecture and alignment data."""
    
    def decide(
        self, 
        architecture_result: Dict[str, Any], 
        alignment_result: Dict[str, Any]
    ) -> TrimConsensus:
        """Derive consensus trim parameters."""
        
        # Extract Architecture Recommendations
        arch_recs = architecture_result.get("trim_recommendations", {})
        arch_5p = arch_recs.get("recommended_5prime_trim", 0)
        arch_3p = arch_recs.get("recommended_3prime_trim", 0)
        arch_adapter = arch_recs.get("three_prime_adapter")
        arch_match = architecture_result.get("architecture_match")
        
        # Extract Alignment Recommendations
        align_recs = alignment_result.get("trim_recommendations", {})
        align_5p = align_recs.get("recommended_5prime_trim", 0)
        align_3p = align_recs.get("recommended_3prime_trim", 0)
        align_consensus = align_recs.get("consensus_level", 0.0)
        global_align_pattern = align_recs.get("global_pattern_detected", False)

        # 1. 5' Trim Decision
        final_5p = 0
        method_5p = "unknown"
        
        if arch_5p == align_5p:
            # perfect agreement
            final_5p = arch_5p
            method_5p = "agreement"
        elif arch_match and arch_5p > 0:
            # Trust known architecture over alignment noise if explicit match
            final_5p = arch_5p
            method_5p = "architecture_dominant"
        elif global_align_pattern and align_consensus > 0.8:
            # Strong alignment signal overrides generic architecture
            final_5p = align_5p
            method_5p = "alignment_dominant"
        else:
            # Fallback to safer option (usually architecture if detected, or 0)
            final_5p = arch_5p if arch_match else align_5p
            method_5p = "fallback"

        # 2. 3' Trim Decision (Adapter vs Fixed)
        final_3p = 0
        final_adapter = arch_adapter
        
        # Comparison Logic
        details = {
            "sources": {
                "architecture": {"5p": arch_5p, "3p": arch_3p, "match": arch_match},
                "alignment": {"5p": align_5p, "3p": align_3p, "consensus": align_consensus}
            },
            "decisions": {
                "5p_logic": method_5p
            }
        }
        
        return TrimConsensus(
            trim_5p=final_5p,
            trim_3p=final_3p, # Usually 0 if adapter handled
            adapter_sequence=final_adapter,
            confidence=0.9 if method_5p == "agreement" else 0.7,
            method="combined_consensus",
            details=details
        )
