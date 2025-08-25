"""RPF-specific check implementations for getRPF.

This module provides specialized checks for Ribosome Protected Fragment (RPF) data.
Each check analyzes specific aspects of the sequence data to determine if it
matches expectations for clean RPF data.
"""

import logging
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, Any

from .processors.check import CleanlinessResults

logger = logging.getLogger(__name__)


class Status(Enum):
    """Status of a check result."""

    PASS = "PASS"
    WARN = "WARNING"
    FAIL = "FAIL"


@dataclass
class CheckResult:
    """Result from a single check.

    Attributes:
        status: Overall status (PASS/WARN/FAIL)
        message: Explanation of the result
        details: Additional numerical/statistical details
    """

    status: Status
    message: str
    details: Optional[Dict] = None


class BaseCheck:
    """Base class for all RPF-specific checks."""

    def check(self, results: CleanlinessResults) -> CheckResult:
        """Run check on CleanlinessResults.

        Args:
            results: CleanlinessResults object containing sequence statistics

        Returns:
            CheckResult containing status and details
        """
        raise NotImplementedError("Subclasses must implement check()")


class LengthDistributionCheck(BaseCheck):
    """Check if length distribution matches RPF expectations."""

    def __init__(
        self,
        min_length: int = 26,
        max_length: int = 35,
        warn_threshold: float = 0.5,
        fail_threshold: float = 0.3,
    ):
        """Initialize length distribution check.

        Args:
            min_length: Minimum expected RPF length
            max_length: Maximum expected RPF length
            warn_threshold: Warn if fraction in range below this
            fail_threshold: Fail if fraction in range below this
        """
        self.min_length = min_length
        self.max_length = max_length
        self.warn_threshold = warn_threshold
        self.fail_threshold = fail_threshold

    def check(self, results: CleanlinessResults) -> CheckResult:
        """Check if read length distribution matches RPF expectations."""
        dist = results.length_distribution
        total_reads = sum(dist.values())

        # Calculate reads in RPF range
        in_range = sum(
            count
            for length, count in dist.items()
            if self.min_length <= length <= self.max_length
        )
        fraction_in_range = in_range / total_reads if total_reads > 0 else 0

        # Find mode
        mode_length = max(dist.items(), key=lambda x: x[1])[0]
        mode_fraction = dist[mode_length] / total_reads

        # Determine status
        if fraction_in_range < self.fail_threshold:
            status = Status.FAIL
            message = (
                f"Only {fraction_in_range:.1%} of reads in RPF range "
                f"({self.min_length}-{self.max_length})"
            )
        elif fraction_in_range < self.warn_threshold:
            status = Status.WARN
            message = (
                f"Low fraction ({fraction_in_range:.1%}) of reads in RPF range "
                f"({self.min_length}-{self.max_length})"
            )
        elif not (self.min_length <= mode_length <= self.max_length):
            status = Status.WARN
            message = (
                f"Mode length ({mode_length}) outside expected RPF range "
                f"({self.min_length}-{self.max_length})"
            )
        else:
            status = Status.PASS
            message = "Read length distribution matches RPF expectations"

        return CheckResult(
            status=status,
            message=message,
            details={
                "fraction_in_range": fraction_in_range,
                "mode_length": mode_length,
                "mode_fraction": mode_fraction,
                "total_reads": total_reads,
                "min_length": self.min_length,
                "max_length": self.max_length,
            },
        )


class BaseCompositionCheck(BaseCheck):
    """Check for concerning patterns in base composition."""

    def __init__(self, max_base_freq: float = 0.85):
        """Initialize base composition check.

        Args:
            max_base_freq: Maximum allowable frequency for any base at a position
        """
        self.max_base_freq = max_base_freq

    def check(self, results: CleanlinessResults) -> CheckResult:
        """Check for positions with extreme base bias."""
        problems = []
        max_freq = 0.0
        worst_pos = None
        worst_base = None

        # Check each position
        positions = len(next(iter(results.nucleotide_frequencies.values())))
        for pos in range(positions):
            pos_freqs = {
                base: freqs[pos]
                for base, freqs in results.nucleotide_frequencies.items()
            }
            pos_max_freq = max(pos_freqs.values())
            if pos_max_freq > self.max_base_freq:
                max_base = max(pos_freqs.items(), key=lambda x: x[1])[0]
                problems.append(f"Position {pos}: {max_base}={pos_max_freq:.1%}")
                if pos_max_freq > max_freq:
                    max_freq = pos_max_freq
                    worst_pos = pos
                    worst_base = max_base

        # Determine status
        if problems:
            status = Status.FAIL
            message = (
                f"Found {len(problems)} positions with extreme base bias "
                f"(>{self.max_base_freq:.1%})"
            )
        else:
            status = Status.PASS
            message = "Base composition looks normal"

        return CheckResult(
            status=status,
            message=message,
            details={
                "problem_positions": problems,
                "max_frequency": max_freq,
                "worst_position": worst_pos,
                "worst_base": worst_base,
            },
        )


class GCContentCheck(BaseCheck):
    """Check if GC content is within expected range."""

    def __init__(
        self,
        min_gc: float = 0.35,
        max_gc: float = 0.65,
    ):
        """Initialize GC content check.

        Args:
            min_gc: Minimum expected GC content (fraction)
            max_gc: Maximum expected GC content (fraction)
        """
        self.min_gc = min_gc
        self.max_gc = max_gc

    def check(self, results: CleanlinessResults) -> CheckResult:
        """Check if GC content is within expected range."""
        if results.gc_content is None:
            return CheckResult(
                status=Status.WARN,
                message="GC content not available",
                details={"gc_content": None},
            )

        gc_fraction = results.gc_content / 100  # Convert from percentage

        if gc_fraction < self.min_gc:
            status = Status.WARN
            message = f"GC content ({gc_fraction:.1%}) below expected minimum ({self.min_gc:.1%})"
        elif gc_fraction > self.max_gc:
            status = Status.WARN
            message = f"GC content ({gc_fraction:.1%}) above expected maximum ({self.max_gc:.1%})"
        else:
            status = Status.PASS
            message = "GC content within expected range"

        return CheckResult(
            status=status, message=message, details={"gc_content": gc_fraction}
        )


class InformationContentCheck(BaseCheck):
    """Check for consistent information content across read positions.
    
    This is the primary check for detecting adapter contamination and repetitive
    sequences by measuring Shannon entropy at each position.
    """
    
    def __init__(self, min_entropy: float = 1.0, ignore_end_positions: int = 2):
        """Initialize information content check.
        
        Args:
            min_entropy: Minimum Shannon entropy threshold per position
            ignore_end_positions: Number of positions to ignore at each end
        """
        self.min_entropy = min_entropy
        self.ignore_end_positions = ignore_end_positions
    
    def check(self, results: CleanlinessResults) -> CheckResult:
        """Check Shannon entropy at each position for uniform complexity."""
        entropies = []
        low_entropy_positions = []
        
        positions = len(next(iter(results.nucleotide_frequencies.values())))
        
        # Only check middle positions, ignore ends
        start_pos = self.ignore_end_positions
        end_pos = positions - self.ignore_end_positions
        
        for pos in range(start_pos, end_pos):
            # Get frequencies at this position
            freqs = [results.nucleotide_frequencies[nt][pos] for nt in "ACGT"]
            # Filter out zero frequencies for entropy calculation
            freqs = [f for f in freqs if f > 0]
            
            if len(freqs) <= 1:
                entropy = 0.0  # No diversity
            else:
                # Calculate Shannon entropy
                entropy = -sum(f * math.log2(f) for f in freqs)
            
            entropies.append(entropy)
            
            if entropy < self.min_entropy:
                low_entropy_positions.append((pos, entropy))
        
        mean_entropy = sum(entropies) / len(entropies) if entropies else 0
        min_entropy = min(entropies) if entropies else 0
        
        # Determine status
        if low_entropy_positions:
            status = Status.FAIL
            message = f"Found {len(low_entropy_positions)} positions with low information content (entropy < {self.min_entropy})"
        else:
            status = Status.PASS
            message = "Information content uniform across positions"
        
        return CheckResult(
            status=status,
            message=message, 
            details={
                "mean_entropy": mean_entropy,
                "min_entropy": min_entropy,
                "low_entropy_positions": low_entropy_positions,
                "entropy_threshold": self.min_entropy
            }
        )


class EndBiasCheck(BaseCheck):
    """Check for extreme nucleotide bias at 5' and 3' ends of reads.
    
    Note: This detects single-nucleotide bias (e.g. poly-A tails, primer sites).
    For adapter contamination (repetitive sequences), use InformationContentCheck.
    """
    
    def __init__(self, end_positions: int = 3, max_bias: float = 0.7):
        """Initialize end bias check.
        
        Args:
            end_positions: Number of positions to check at each end
            max_bias: Maximum allowable nucleotide frequency at ends
        """
        self.end_positions = end_positions
        self.max_bias = max_bias
    
    def check(self, results: CleanlinessResults) -> CheckResult:
        """Check for excessive nucleotide bias at read ends.
        
        For 3' end analysis, we use reversed sequences to align all reads at their 3' ends,
        making adapter contamination detectable as nucleotide bias.
        """
        problems = []
        max_bias_found = 0.0
        worst_position = None
        
        positions = len(next(iter(results.nucleotide_frequencies.values())))
        
        # Check 5' end positions (straightforward - all reads start at pos 0)
        for pos in range(min(self.end_positions, positions)):
            pos_freqs = {nt: results.nucleotide_frequencies[nt][pos] for nt in "ACGT"}
            max_freq = max(pos_freqs.values())
            
            if max_freq > self.max_bias:
                dominant_nt = max(pos_freqs.items(), key=lambda x: x[1])[0]
                problems.append(f"5' position {pos}: {dominant_nt}={max_freq:.1%}")
                if max_freq > max_bias_found:
                    max_bias_found = max_freq
                    worst_position = f"5'_{pos}"
        
        # Check 3' end positions using REVERSED sequences
        if hasattr(results, 'reversed_nucleotide_frequencies') and results.reversed_nucleotide_frequencies:
            reversed_positions = len(next(iter(results.reversed_nucleotide_frequencies.values())))
            
            for pos in range(min(self.end_positions, reversed_positions)):
                pos_freqs = {nt: results.reversed_nucleotide_frequencies[nt][pos] for nt in "ACGT"}
                max_freq = max(pos_freqs.values())
                
                if max_freq > self.max_bias:
                    dominant_nt = max(pos_freqs.items(), key=lambda x: x[1])[0]
                    problems.append(f"3' position {pos} (reversed): {dominant_nt}={max_freq:.1%}")
                    if max_freq > max_bias_found:
                        max_bias_found = max_freq
                        worst_position = f"3'_reversed_{pos}"
        else:
            # Fallback to old method if reversed frequencies not available
            for i in range(min(self.end_positions, positions)):
                pos = positions - 1 - i  # Count from end
                pos_freqs = {nt: results.nucleotide_frequencies[nt][pos] for nt in "ACGT"}
                max_freq = max(pos_freqs.values())
                
                if max_freq > self.max_bias:
                    dominant_nt = max(pos_freqs.items(), key=lambda x: x[1])[0]
                    problems.append(f"3' position -{i+1}: {dominant_nt}={max_freq:.1%}")
                    if max_freq > max_bias_found:
                        max_bias_found = max_freq
                        worst_position = f"3'_{i+1}"
        
        # Determine status
        if problems:
            status = Status.FAIL
            message = f"Found end bias at {len(problems)} positions (>{self.max_bias:.1%} frequency)"
        else:
            status = Status.PASS
            message = "No excessive nucleotide bias at read ends"
        
        return CheckResult(
            status=status,
            message=message,
            details={
                "biased_positions": problems,
                "max_bias_found": max_bias_found,
                "worst_position": worst_position,
                "bias_threshold": self.max_bias,
                "uses_reversed_analysis": hasattr(results, 'reversed_nucleotide_frequencies')
            }
        )


class SoftClippingCheck(BaseCheck):
    """Check for excessive soft-clipping in alignments."""
    
    def __init__(self, max_clip_rate: float = 0.1, max_mean_clips: float = 1.0):
        """Initialize soft-clipping check.
        
        Args:
            max_clip_rate: Maximum fraction of reads with soft clips
            max_mean_clips: Maximum mean soft clips per read
        """
        self.max_clip_rate = max_clip_rate
        self.max_mean_clips = max_mean_clips
    
    def check(self, alignment_stats: Dict[str, Any]) -> CheckResult:
        """Check soft-clipping statistics from alignment."""
        if "total_reads_analyzed" not in alignment_stats:
            return CheckResult(
                status=Status.WARN,
                message="Soft-clipping statistics not available",
                details={}
            )
        
        total_reads = alignment_stats["total_reads_analyzed"]
        reads_with_clips = alignment_stats["reads_with_soft_clips"]
        mean_5prime = alignment_stats["mean_5prime_clips"]
        mean_3prime = alignment_stats["mean_3prime_clips"]
        
        if total_reads == 0:
            clip_rate = 0
        else:
            clip_rate = reads_with_clips / total_reads
        
        mean_total_clips = mean_5prime + mean_3prime
        
        problems = []
        if clip_rate > self.max_clip_rate:
            problems.append(f"High soft-clipping rate: {clip_rate:.1%}")
        
        if mean_total_clips > self.max_mean_clips:
            problems.append(f"High mean soft clips: {mean_total_clips:.1f}")
        
        if problems:
            status = Status.FAIL
            message = "; ".join(problems)
        else:
            status = Status.PASS
            message = "Soft-clipping within acceptable limits"
        
        return CheckResult(
            status=status,
            message=message,
            details={
                "clip_rate": clip_rate,
                "mean_total_clips": mean_total_clips,
                "reads_with_clips": reads_with_clips,
                "total_reads": total_reads
            }
        )


def categorize_failures(results: Dict[str, CheckResult]) -> Dict[str, str]:
    """Categorize sample by failure type for seqspec batch processing.
    
    Args:
        results: Dictionary mapping check names to results
        
    Returns:
        Dictionary with failure categories
    """
    failure_categories = []
    
    for check_name, result in results.items():
        if result.status == Status.FAIL:
            if "length" in check_name.lower():
                failure_categories.append("length_distribution")
            elif "end" in check_name.lower() or "bias" in check_name.lower():
                failure_categories.append("end_bias")
            elif "information" in check_name.lower() or "entropy" in check_name.lower():
                failure_categories.append("low_complexity")
            elif "clipping" in check_name.lower():
                failure_categories.append("soft_clipping")
            elif "composition" in check_name.lower():
                failure_categories.append("base_composition")
            else:
                failure_categories.append("other")
    
    return {
        "failure_categories": failure_categories,
        "primary_failure": failure_categories[0] if failure_categories else None,
        "is_clean": len(failure_categories) == 0
    }


def write_check_report(results: Dict[str, CheckResult], output: Path) -> None:
    """Write check results to a file.

    Args:
        results: Dictionary mapping check names to results
        output: Path where to write the report
    """
    with open(output, "w") as f:
        # Write summary
        f.write("=== RPF Data Check Results ===\n\n")
        
        # Add cleanliness status
        categories = categorize_failures(results)
        f.write(f"Sample Status: {'CLEAN' if categories['is_clean'] else 'NEEDS_SEQSPEC'}\n")
        if not categories['is_clean']:
            f.write(f"Primary Failure Type: {categories['primary_failure']}\n")
            f.write(f"All Failure Types: {', '.join(categories['failure_categories'])}\n")
        f.write("\n")
        
        f.write("Check Summary:\n")
        for check_name, result in results.items():
            f.write(f"{check_name:30} [{result.status.value}]\n")

        # Write details
        f.write("\nDetailed Results:\n")
        for check_name, result in results.items():
            f.write(f"\n{check_name}:\n")
            f.write(f"Status: {result.status.value}\n")
            f.write(f"Message: {result.message}\n")

            if result.details:
                f.write("Details:\n")
                for key, value in result.details.items():
                    if isinstance(value, list):
                        f.write(f"  {key}:\n")
                        for item in value:
                            f.write(f"    - {item}\n")
                    else:
                        f.write(f"  {key}: {value}\n")


def run_all_cleanliness_checks(
    sequence_results: CleanlinessResults,
    alignment_stats: Optional[Dict[str, Any]] = None
) -> Dict[str, CheckResult]:
    """Run all cleanliness checks and return results.
    
    Args:
        sequence_results: Results from sequence analysis
        alignment_stats: Optional alignment statistics
        
    Returns:
        Dictionary mapping check names to results
    """
    checks = {
        "length_distribution": LengthDistributionCheck(),
        "information_content": InformationContentCheck(),
        "end_bias": EndBiasCheck(),
        "base_composition": BaseCompositionCheck(),
        "gc_content": GCContentCheck()
    }
    
    results = {}
    
    # Run sequence-based checks
    for check_name, check in checks.items():
        results[check_name] = check.check(sequence_results)
    
    # Run alignment-based checks if available
    if alignment_stats:
        soft_clip_check = SoftClippingCheck()
        results["soft_clipping"] = soft_clip_check.check(alignment_stats)
    
    return results
