"""Core processing modules for getRPF."""

from .adapter import AdapterDetector, analyze_adapter_file
from .check import CleanlinessChecker, CleanlinessResults, analyze_file

__all__ = [
    "AdapterDetector",
    "CleanlinessChecker",
    "CleanlinessResults",
    "analyze_adapter_file",
    "analyze_file",
]
