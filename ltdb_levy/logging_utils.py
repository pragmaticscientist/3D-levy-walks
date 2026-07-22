"""Structured logging setup.

The pipeline makes a lot of decisions that must remain auditable after the fact
(which files were excluded, which transform was applied to which coordinates,
which conditions failed the pool minimums). Those go through this logger so they
land in both the console and a per-run log file.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

_LOG_FORMAT = "%(asctime)s %(levelname)-7s %(name)-28s %(message)s"
_DATE_FORMAT = "%H:%M:%S"


def configure_logging(
    level: str = "INFO",
    log_file: Optional[Path] = None,
) -> None:
    """Configure the root logger for a pipeline run.

    Args:
        level: Logging level name, e.g. ``"INFO"`` or ``"DEBUG"``.
        log_file: If given, also write the full log to this path. The parent
            directory is created if needed.
    """
    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Replace handlers so repeated CLI invocations in one process do not
    # accumulate duplicate output.
    for handler in list(root.handlers):
        root.removeHandler(handler)

    formatter = logging.Formatter(_LOG_FORMAT, datefmt=_DATE_FORMAT)

    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    root.addHandler(stream_handler)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    """Return the module-scoped logger for ``name``."""
    return logging.getLogger(name)
