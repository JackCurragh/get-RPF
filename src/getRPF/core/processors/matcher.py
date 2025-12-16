"""Strict architecture matching logic.

This module validates read architectures against observed signal statistics
using deterministic boolean logic constraints.
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Tuple
from .signals import SignalStats
from .types import ReadArchitecture

logger = logging.getLogger(__name__)


@dataclass
class MatchResult:
    """Result of an architecture match attempt."""
    architecture: ReadArchitecture
    is_match: bool
    confidence: float
    reasons: List[str]  # Explanations for acceptance/rejection


class ArchitectureMatcher:
    """Matches partial read signals to known architectures using strict rules."""

    def __init__(self, confidence_threshold: float = 0.8):
        self.threshold = confidence_threshold

    def match(self, stats: SignalStats, architectures: List[ReadArchitecture]) -> Optional[MatchResult]:
        """Find the best matching architecture from a list.
        
        Args:
            stats: Signal statistics from the sample
            architectures: List of candidate architectures
            
        Returns:
            Best matching result or None if no match found
        """
        best_result = None
        best_score = 0.0
        
        for arch in architectures:
            result = self._evaluate_architecture(stats, arch)
            
            # Log the decision trace
            status = "MATCH" if result.is_match else "REJECT"
            logger.debug(f"Architecture {arch.protocol_name}: {status} (Score: {result.confidence:.2f})")
            for reason in result.reasons:
                logger.debug(f"  - {reason}")
            
            if result.is_match and result.confidence > best_score:
                best_score = result.confidence
                best_result = result
        
        return best_result

    def _evaluate_architecture(self, stats: SignalStats, arch: ReadArchitecture) -> MatchResult:
        """Evaluate a single architecture against the signals."""
        reasons = []
        score = 0.0
        
        # 1. Adapter Constraint (Critical)
        # If an architecture specifies an adapter, we MUST see it.
        # However, we don't have raw reads here, just stats.
        # Wait - strict matching usually requires checking specific sequences in potential adapter regions.
        # The SignalStats gives us entropy/composition. 
        # Low entropy at specific positions suggests adapter.
        # But to be robust, we probably need the raw reads for the adapter check 
        # or we accept that the SignalProcessor 'entropy' drop is the proxy.
        # BETTER: The Matcher should ideally have access to 'Adapter Presence' boolean 
        # which might be calculated separately?
        # OR: We keep the fuzzy adapter check as a pre-filter?
        # NO, the plan said "Strict Matcher".
        # Let's assume we pass in 'adapter_presence_score' computed externally?
        # Actually, let's implement a targeted adapter check here if we can, 
        # but we only have stats.
        
        # Checking UMI Profile (Entropy Check)
        umi_score = self._check_umi_profile(stats, arch, reasons)
        
        # Checking Length Distribution
        # We don't have the length distribution in SignalStats directly (we have max len),
        # but we can infer RPF length compatibility?
        # Actually, SignalStats needs length distribution to be truly useful here.
        # Let's assume we can get simple length checks from valid_reads passed elsewhere,
        # or we update SignalStats.
        # For this iteration, let's focus on the signals we HAVE.
        
        # Length check (placeholder logic based on SignalStats max_len unfortunately)
        # Real implementation should probably pass length hist in SignalStats.
        
        # Let's implement the logic based on entropy profile matching.
        
        # Score Accumulation
        # If UMI check failed (and UMI was expected), it's a hard reject.
        if umi_score == 0.0 and arch.umi_positions:
            return MatchResult(arch, False, 0.0, reasons)
            
        score += umi_score * 0.4
        
        # Adapter check needs to be handled. 
        # Ideally, we'd check if specific adapter k-mers are overrepresented in the composition stats.
        # e.g. if adapter is "AGATCGG", check if pos X has 'A', pos X+1 has 'G' etc.
        # This is actually very robust!
        adapter_score = self._check_adapter_composition(stats, arch, reasons)
        if adapter_score == 0.0 and arch.adapter_sequences:
             return MatchResult(arch, False, 0.0, reasons)
             
        score += adapter_score * 0.6
        
        is_match = score >= self.threshold
        
        return MatchResult(arch, is_match, score, reasons)

    def _check_umi_profile(self, stats: SignalStats, arch: ReadArchitecture, reasons: List[str]) -> float:
        """Verify UMI entropy profile."""
        if not arch.umi_positions:
            reasons.append("No UMI expected (Pass)")
            return 1.0
            
        for start, end in arch.umi_positions:
            # Check entropy in this region
            # UMI region should have high entropy (> 1.5 usually)
            # Need to be careful about 5' vs 3' coordinates. 
            # arch.umi_positions are usually 5' relative.
            
            region_entropy = stats.entropy_5p[start:end]
            if not region_entropy:
                 reasons.append(f"UMI region {start}-{end} out of bounds")
                 return 0.0
                 
            avg_entropy = sum(region_entropy) / len(region_entropy)
            if avg_entropy < 1.0:
                reasons.append(f"UMI region {start}-{end} has low entropy ({avg_entropy:.2f} < 1.0)")
                return 0.0
            else:
                reasons.append(f"UMI region {start}-{end} has high entropy ({avg_entropy:.2f})")
                
        return 1.0

    def _check_adapter_composition(self, stats: SignalStats, arch: ReadArchitecture, reasons: List[str]) -> float:
        """Check if consensus composition matches expected adapter."""
        if not arch.adapter_sequences:
            return 1.0
            
        # We can try to see if the adapter sequence 'emerges' in the 3' composition stats.
        # 3' adapter means the read ends with it (or we trimmed it?).
        # If we trimmed it, it won't be there. 
        # If we are detecting it in raw reads, it should be at the 3' end if the read is short,
        # or present in the sequence.
        
        # For "Ingolia 2009": 3' linker "CTGTAG..."
        # In `composition_3p`, index 0 is the LAST base.
        # If adapter is at the end, `composition_3p[0]` should resemble the last base of adapter check?
        # Actually, adapter detection usually happens BEFORE this matching logic in the old code.
        # But here we want strict evidence.
        
        # Strategy: Check if any of the architectures adapters matches the consensus sequence
        # derived from the composition stats, EITHER at 5' (rare) or 3' (common).
        
        best_adapter_score = 0.0
        
        for adapter in arch.adapter_sequences:
            # Check 3' end match (most common)
            # Reverse adapter to match 3' stats (index 0 = last base)
            rev_adapter = adapter[::-1]
            match_len = min(len(rev_adapter), len(stats.composition_3p))
            matches = 0
            
            for i in range(match_len):
                expected_base = rev_adapter[i]
                # Check frequency of this base at pos i from 3' end
                if stats.composition_3p[i].get(expected_base, 0) > 0.4: # Consensus threshold
                    matches += 1
            
            # If we see > 50% of the adapter bases in the consensus, that's a strong signal
            # for un-trimmed data.
            # However, if data is already trimmed, this will fail.
            # We need to handle that ambiguity?
            # Assuming Raw Data for architecture detection.
            
            score = matches / match_len if match_len > 0 else 0
            if score > best_adapter_score:
                best_adapter_score = score
                
        if best_adapter_score > 0.7:
             reasons.append(f"Strong 3' adapter signal found (Score: {best_adapter_score:.2f})")
             return 1.0
        elif best_adapter_score > 0.4:
             reasons.append(f"Weak 3' adapter signal found (Score: {best_adapter_score:.2f})")
             return 0.5
        else:
             # It's possible the adapter is variable position (internal), not strictly at 3' end.
             # Or it's not present (clean data?).
             # For strict matching, if we DON'T see it, we might reject if we expect raw data.
             reasons.append(f"Adapter footprint not found at 3' end (Score: {best_adapter_score:.2f})")
             return 0.0
