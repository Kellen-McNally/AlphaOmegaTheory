"""
Centralized logging configuration for the αΩ Framework.

This module provides a consistent logging interface across all modules.
It supports:
- Multiple log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- Console and file output
- Formatted output with timestamps and module names
- Easy integration into existing code

Usage in any module:
    from utils.logging_config import get_logger
    logger = get_logger(__name__)

    logger.debug("Detailed information for debugging")
    logger.info("General informational messages")
    logger.warning("Warning messages")
    logger.error("Error messages")
    logger.critical("Critical errors")
"""

import logging
import sys
from pathlib import Path
from typing import Optional


# Default log level for the framework
DEFAULT_LOG_LEVEL = logging.INFO

# Log format: timestamp - module - level - message
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
DATE_FORMAT = '%Y-%m-%d %H:%M:%S'

# Global flag to track if logging has been configured
_logging_configured = False


def configure_logging(
    level: int = DEFAULT_LOG_LEVEL,
    log_file: Optional[str] = None,
    console: bool = True,
    format_string: Optional[str] = None,
    date_format: Optional[str] = None
) -> None:
    """
    Configure the root logger for the αΩ Framework.

    This should typically be called once at the start of a script or application.
    If not called explicitly, get_logger() will use default console logging.

    Args:
        level: Logging level (logging.DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file. If provided, logs will be written to file.
        console: If True, output to console (stderr). Default True.
        format_string: Optional custom format string. Uses LOG_FORMAT if not provided.
        date_format: Optional custom date format. Uses DATE_FORMAT if not provided.

    Example:
        # Console only, INFO level
        configure_logging(level=logging.INFO)

        # Console + file, DEBUG level
        configure_logging(level=logging.DEBUG, log_file='aomega.log')

        # File only, WARNING level
        configure_logging(level=logging.WARNING, log_file='warnings.log', console=False)
    """
    global _logging_configured

    # Use default formats if not provided
    if format_string is None:
        format_string = LOG_FORMAT
    if date_format is None:
        date_format = DATE_FORMAT

    # Create formatter
    formatter = logging.Formatter(format_string, datefmt=date_format)

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove existing handlers to avoid duplicates
    root_logger.handlers.clear()

    # Add console handler if requested
    if console:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)

    # Add file handler if requested
    if log_file:
        # Create log directory if needed
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    _logging_configured = True


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for the specified module.

    This is the main function to use in modules throughout the framework.
    If configure_logging() hasn't been called, this will set up basic console logging.

    Args:
        name: Module name (typically __name__)

    Returns:
        Logger instance configured for the framework

    Example:
        logger = get_logger(__name__)
        logger.info("Starting calculation...")
    """
    global _logging_configured

    # If logging not configured, set up basic console logging
    if not _logging_configured:
        configure_logging(level=DEFAULT_LOG_LEVEL, console=True)

    return logging.getLogger(name)


def set_log_level(level: int) -> None:
    """
    Change the log level for all loggers.

    Args:
        level: New logging level (logging.DEBUG, INFO, WARNING, ERROR, CRITICAL)

    Example:
        set_log_level(logging.DEBUG)  # Enable debug output
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Update all handlers
    for handler in root_logger.handlers:
        handler.setLevel(level)


def disable_logging() -> None:
    """
    Disable all logging output.

    Useful for tests or when you want silent operation.
    """
    logging.disable(logging.CRITICAL)


def enable_logging() -> None:
    """
    Re-enable logging after disable_logging() was called.
    """
    logging.disable(logging.NOTSET)


# Convenience functions for common logging configurations

def configure_for_cli(verbose: bool = False) -> None:
    """
    Configure logging for CLI/command-line scripts.

    Args:
        verbose: If True, use DEBUG level. Otherwise INFO level.
    """
    level = logging.DEBUG if verbose else logging.INFO
    configure_logging(level=level, console=True)


def configure_for_library() -> None:
    """
    Configure logging for use as a library.

    Sets WARNING level by default so library doesn't spam user's logs.
    Users can call configure_logging() to override.
    """
    configure_logging(level=logging.WARNING, console=True)


def configure_for_tests(log_file: Optional[str] = None) -> None:
    """
    Configure logging for test runs.

    Args:
        log_file: Optional file to capture test logs
    """
    # Tests typically want to see warnings and errors
    configure_logging(level=logging.WARNING, console=True, log_file=log_file)


# Example usage and self-test
if __name__ == '__main__':
    # Configure logging
    configure_logging(level=logging.DEBUG, console=True)

    # Get logger
    logger = get_logger(__name__)

    # Test all levels
    logger.debug("This is a debug message")
    logger.info("This is an info message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")

    print("\n--- Testing log levels ---")

    # Change to INFO level
    set_log_level(logging.INFO)
    logger.debug("This debug message should NOT appear")
    logger.info("This info message SHOULD appear")

    print("\n--- Testing file output ---")

    # Configure with file output
    configure_logging(level=logging.INFO, log_file='test_aomega.log', console=True)
    logger = get_logger(__name__)
    logger.info("This message goes to both console and file")

    print("\n[OK] Logging configuration test complete")
    print("Check test_aomega.log for file output")
