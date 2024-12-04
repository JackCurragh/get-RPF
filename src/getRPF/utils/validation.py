"""Validation utilities for getRPF.

This module provides validation functions for various inputs including:
    - File formats
    - Sequence data
    - Configuration parameters
"""

from typing import Set


def validate_input_format(format: str) -> bool:
    """Validate the input format string.

    Args:
        format: Format string to validate

    Returns:
        bool: True if format is valid

    Raises:
        ValueError: If format is not supported

    Examples:
        >>> validate_input_format('fastq')
        True
        >>> validate_input_format('invalid')
        Raises ValueError
    """
    valid_formats: Set[str] = {"fastq", "fasta", "collapsed"}
    if format.lower() not in valid_formats:
        raise ValueError(
            f"Unsupported format: {format}. " f"Must be one of {sorted(valid_formats)}"
        )
    return True


def validate_adapter_sequence(adapter: str) -> bool:
    """Validate the adapter sequence.

    Args:
        adapter: Adapter sequence to validate

    Returns:
        bool: True if adapter sequence is valid

    Raises:
        ValueError: If adapter sequence contains invalid characters

    Examples:
        >>> validate_adapter_sequence('ATCG')
        True
        >>> validate_adapter_sequence('INVALID')
        Raises ValueError
    """
    valid_bases: Set[str] = {"A", "T", "C", "G", "N"}
    if not adapter:
        raise ValueError("Adapter sequence cannot be empty")
    if not all(base.upper() in valid_bases for base in adapter):
        raise ValueError(
            f"Invalid adapter sequence: {adapter}. "
            f"Must contain only {sorted(valid_bases)}"
        )
    return True
