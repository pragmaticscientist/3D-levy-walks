"""Exception hierarchy for the ltdb_levy pipeline.

Every failure mode that a user can plausibly trigger gets its own class so the
CLI can report something actionable instead of a bare traceback.
"""

from __future__ import annotations


class LtdbLevyError(Exception):
    """Base class for all errors raised by this package."""


class ConfigError(LtdbLevyError):
    """The YAML configuration is missing, malformed, or internally inconsistent."""


class ParseError(LtdbLevyError):
    """An LTDB ground-truth file could not be parsed."""

    def __init__(self, path: str, line_number: int, message: str) -> None:
        self.path = path
        self.line_number = line_number
        super().__init__("{}:{}: {}".format(path, line_number, message))


class MetadataError(LtdbLevyError):
    """The biological metadata table is missing, incomplete, or unverified."""


class InsufficientDataError(LtdbLevyError):
    """A condition does not meet the configured minimum sample requirements.

    This is raised only where continuing would silently produce a meaningless
    result. Ordinary shortfalls are recorded in an audit table instead.
    """


class FitError(LtdbLevyError):
    """A maximum-likelihood fit failed to converge or was given degenerate input."""


class GeometryError(LtdbLevyError):
    """A target geometry is impossible for the requested projected area."""
