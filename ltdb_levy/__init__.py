"""Fitting Levy-walk surrogates to the Leukocyte Tracking Database (LTDB).

The package parses LTDB ground-truth trajectories, segments them into directed
runs, fits the relocation-length exponent of the simulator's mixture step law,
builds nonparametric relocation-length pools, and replays a ladder of surrogate
trajectory classes against controlled 3D targets.

The power-law hypothesis is tested, never assumed.
"""

from __future__ import annotations

__version__ = "0.1.0"

# Versions of the on-disk artifact schemas. Bump when a column set changes so
# stale intermediate files are detected rather than silently misread.
SCHEMA_VERSIONS = {
    "observations": 1,
    "displacements": 1,
    "runs": 1,
    "pool_audit": 1,
    "fits": 2,
    "replay_trials": 1,
}

__all__ = ["__version__", "SCHEMA_VERSIONS"]
