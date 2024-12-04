"""File handling utilities for getRPF.

This module provides utilities for:
    - File validation and checking
    - File path manipulation
    - Common file operations
"""

import bz2
import gzip
import os
from pathlib import Path
from typing import Union


def check_file_readability(file_path: Union[str, Path]) -> bool:
    """Check if a file exists and is readable.

    Args:
        file_path: Path to file to check

    Returns:
        bool: True if file is readable

    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file cannot be read
        TypeError: If file_path is not str or Path

    Examples:
        >>> check_file_readability('existing_file.txt')
        True
        >>> check_file_readability('nonexistent.txt')
        Raises FileNotFoundError
    """
    if isinstance(file_path, str):
        file_path = Path(file_path)
    elif not isinstance(file_path, Path):
        raise TypeError(f"file_path must be str or Path, not {type(file_path)}")

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if not os.access(file_path, os.R_OK):
        raise PermissionError(f"Cannot read file: {file_path}")
    return True


def get_file_opener(filepath: Path):
    """Determine the appropriate file opener based on file extension."""
    suffix = filepath.suffix.lower()
    if suffix == ".gz":
        return gzip.open
    elif suffix == ".bz2":
        return bz2.open
    return open
