"""Logging configuration for getRPF.

This module provides logging utilities including:
    - Logging setup and configuration
    - Custom formatters and handlers
    - Logging context managers
"""

import logging
from pathlib import Path
from typing import List, Optional


def setup_logging(
    log_level: str = "INFO",
    log_file: Optional[Path] = None,
    module_name: str = "getRPF",
) -> None:
    """Configure logging for the module.

    Sets up logging with consistent formatting and optional file output.

    Args:
        log_level: Logging level (default: "INFO")
        log_file: Optional path to log file
        module_name: Name to use for logger (default: "getRPF")

    Raises:
        ValueError: If log_level is invalid
        PermissionError: If log_file cannot be written

    Example:
        >>> setup_logging(log_level="DEBUG", log_file=Path("app.log"))

    Notes:
        - Log format includes timestamp, logger name, level, and message
        - Stream handler (stdout) is always added
        - File handler is optional
    """
    # Validate log level
    try:
        numeric_level = getattr(logging, log_level.upper())
    except AttributeError:
        raise ValueError(f"Invalid log level: {log_level}")

    # Configure handlers
    handlers: List[logging.Handler] = [logging.StreamHandler()]
    if log_file is not None:
        try:
            handlers.append(logging.FileHandler(log_file))
        except PermissionError:
            raise PermissionError(f"Cannot write to log file: {log_file}")

    # Setup basic configuration
    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )

    # Get logger for module
    logger = logging.getLogger(module_name)
    logger.setLevel(numeric_level)
