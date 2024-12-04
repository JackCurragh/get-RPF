"""RPF-specific check implementations for getRPF.

This module provides specialized checks for Ribosome Protected Fragment (RPF) data.
Each check analyzes specific aspects of the sequence data to determine if it
matches expectations for clean RPF data.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Dict, Optional

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


def write_check_report(results: Dict[str, CheckResult], output: Path) -> None:
    """Write check results to a file.

    Args:
        results: Dictionary mapping check names to results
        output: Path where to write the report
    """
    with open(output, "w") as f:
        # Write summary
        f.write("=== RPF Data Check Results ===\n\n")
        f.write("Summary:\n")
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
