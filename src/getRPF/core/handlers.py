"""Core handlers for getRPF functionality.

This module implements the main processing logic for getRPF commands, serving as an
intermediary layer between the CLI and the core processing modules.

The handlers in this module follow a consistent pattern:
1. Validate inputs and options
2. Initialize appropriate processor objects
3. Process data in a memory-efficient manner
4. Aggregate and format results
5. Clean up temporary resources
"""

import logging
from pathlib import Path
from typing import Optional

from ..core.checkers import (
    BaseCompositionCheck,
    GCContentCheck,
    LengthDistributionCheck,
    write_check_report,
)
from ..core.processors.adapter import AdapterDetector
from ..core.processors.check import CleanlinessChecker
from ..utils.validation import validate_adapter_sequence

# Configure logging
logger = logging.getLogger(__name__)


def handle_cleanliness_check(
    input_file: Path,
    format: str,
    output: Path,
    count_pattern: Optional[str] = None,
    min_quality: Optional[int] = 20,
    threads: Optional[int] = 1,
    rpf_checks: bool = True,
) -> None:
    """Handle the quality check command workflow.

    Args:
        input_file: Path to input sequence file
        format: Format of input file (fastq/fasta/collapsed)
        output: Path where report will be written
        count_pattern: For collapsed format, pattern to extract read counts
        min_quality: Minimum quality score threshold
        threads: Number of processing threads
        rpf_checks: Whether to run RPF-specific checks
    """
    logger.info(f"Starting cleanliness check for {input_file}")

    try:
        # Initialize checker with parameters
        checker = CleanlinessChecker(
            format=format,
            min_quality=min_quality,
            threads=threads,
            count_pattern=count_pattern if format == "collapsed" else None,
        )

        # Process file and get basic results
        results = checker.analyze_file(input_file)

        # Write basic report
        results.write_report(output)

        # Run RPF-specific checks if requested
        if rpf_checks:
            checks = {
                "Length Distribution": LengthDistributionCheck(),
                "Base Composition": BaseCompositionCheck(),
                "GC Content": GCContentCheck(),
            }

            # Run checks
            check_results = {
                name: check.check(results) for name, check in checks.items()
            }

            # Write RPF check report
            rpf_report = output.with_suffix(".rpf_checks.txt")
            write_check_report(check_results, rpf_report)
            logger.info(f"RPF check report written to {rpf_report}")

        logger.info("Cleanliness check completed successfully")

    except Exception as e:
        logger.error(f"Cleanliness check failed: {str(e)}")
        raise RuntimeError(f"Cleanliness check failed: {str(e)}") from e


def handle_adapter_detection(
    input_file: Path,
    format: str,
    output: Path,
    adapter: str,
    min_overlap: Optional[int] = 10,
    max_mismatches: Optional[int] = 1,
    threads: Optional[int] = 1,
    count_pattern: Optional[str] = None,
) -> None:
    """Handle the adapter detection command workflow.

    Args:
        input_file: Path to input sequence file
        format: Format of input file (fastq/fasta/collapsed)
        output: Path where report will be written
        adapter: Sequence of adapter to detect
        min_overlap: Minimum overlap for adapter detection
        max_mismatches: Maximum allowed mismatches
        threads: Number of processing threads
        count_pattern: Pattern to count in header
    """
    logger.info(f"Starting adapter detection for {input_file}")

    try:
        # Validate adapter sequence
        validate_adapter_sequence(adapter)

        # Initialize detector with parameters
        detector = AdapterDetector(
            adapter=adapter,
            format=format,
            min_overlap=min_overlap,
            max_mismatches=max_mismatches,
            threads=threads,
            count_pattern=count_pattern,
        )

        # Process file and get results
        results = detector.analyze_file(input_file)

        # Write report
        results.write_report(output)

        logger.info("Adapter detection completed successfully")

    except Exception as e:
        logger.error(f"Adapter detection failed: {str(e)}")
        raise RuntimeError(f"Adapter detection failed: {str(e)}") from e
