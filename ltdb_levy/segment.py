"""Directed-run segmentation -- the primary definition of an empirical relocation.

Consecutive moving displacements are merged into one *run*. The current run ends
when any of the following occurs:

1. the fragment ends (a frame gap or the end of the track);
2. the step is not "moving" (zero displacement, or speed below threshold);
3. the turning angle relative to the previous moving step exceeds the threshold.

Turning-point convention (deterministic, and tested)
----------------------------------------------------
The displacement whose turning angle *exceeds* the threshold **starts the new
run**; it is never attributed to the old one. Comparisons are strict (``>``), so
an angle exactly equal to the threshold does not split. Where several candidate
split points are adjacent, the lowest step index wins by construction, since the
scan is left to right.

The relocation length is the **end-to-end** norm of the run's net displacement,
not the summed path length. Both are stored, so straightness is available as a
diagnostic.

The first and last run of every fragment are flagged: they are censored by the
cell entering or leaving the field of view, and are excluded from the primary
fit (and included in a sensitivity analysis).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .logging_utils import get_logger

logger = get_logger(__name__)

RUN_COLUMNS: Tuple[str, ...] = (
    "file",
    "cohort",
    "video",
    "population",
    "video_id",
    "track_id",
    "track_uid",
    "fragment_uid",
    "run_id",
    "run_uid",
    "step_start",
    "step_end",
    "n_steps",
    "frame_from",
    "frame_to",
    "run_vector_x",
    "run_vector_y",
    "run_vector_z",
    "run_length_um",
    "path_length_um",
    "straightness",
    "run_duration_s",
    "mean_speed_um_s",
    "is_first",
    "is_last",
    "included",
    "exclusion_reason",
)


@dataclass
class SegmentParams:
    """Segmentation thresholds.

    Attributes:
        turning_angle_degrees: A turn strictly greater than this ends the run.
        speed_threshold: Steps slower than this are not "moving". ``None``
            disables the speed criterion (zero-length steps still break runs).
        minimum_run_frames: Runs made of fewer steps than this are excluded.
        minimum_run_length: Runs shorter than this (end-to-end, um) are excluded.
        exclude_first_and_last_run: Flag the censored boundary runs for
            exclusion from the primary fit.
    """

    turning_angle_degrees: float = 60.0
    speed_threshold: Optional[float] = None
    minimum_run_frames: int = 1
    minimum_run_length: float = 0.0
    exclude_first_and_last_run: bool = True

    @property
    def turning_angle_rad(self) -> float:
        return float(np.radians(self.turning_angle_degrees))


@dataclass
class SegmentAudit:
    """Counts of runs excluded, by reason."""

    n_runs_total: int = 0
    n_runs_included: int = 0
    excluded_first_last: int = 0
    excluded_too_few_frames: int = 0
    excluded_too_short: int = 0
    non_moving_steps: int = 0
    counts: Dict[str, int] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {"metric": "runs_total", "count": self.n_runs_total},
                {"metric": "runs_included", "count": self.n_runs_included},
                {"metric": "excluded_first_or_last_run", "count": self.excluded_first_last},
                {"metric": "excluded_too_few_frames", "count": self.excluded_too_few_frames},
                {"metric": "excluded_too_short", "count": self.excluded_too_short},
                {"metric": "non_moving_steps", "count": self.non_moving_steps},
            ]
        )


def _moving_mask(
    length: np.ndarray, speed: np.ndarray, speed_threshold: Optional[float]
) -> np.ndarray:
    """A step is moving when it has non-zero length and clears the speed bar."""
    moving = length > 0.0
    if speed_threshold is not None:
        with np.errstate(invalid="ignore"):
            moving &= np.nan_to_num(speed, nan=0.0) >= speed_threshold
    return moving


def _split_indices(
    moving: np.ndarray, angles: np.ndarray, angle_threshold: float
) -> List[Tuple[int, int]]:
    """Return inclusive ``(start, end)`` step-index spans for each run.

    Scans left to right over a single fragment. A run is a maximal span of
    consecutive moving steps, further split wherever the turning angle
    strictly exceeds ``angle_threshold``.
    """
    spans: List[Tuple[int, int]] = []
    n = len(moving)
    i = 0
    while i < n:
        if not moving[i]:
            i += 1
            continue
        start = i
        j = i + 1
        while j < n and moving[j]:
            a = angles[j]
            # NaN angle => undefined (a zero-length step intervened). Since the
            # previous step is moving and contiguous here, NaN cannot occur
            # mid-span, but guard anyway rather than splitting on NaN silently.
            if np.isfinite(a) and a > angle_threshold:
                break
            j += 1
        spans.append((start, j - 1))
        i = j
    return spans


def segment_runs(
    displacements: pd.DataFrame,
    params: SegmentParams,
) -> Tuple[pd.DataFrame, SegmentAudit]:
    """Segment frame displacements into directed runs.

    Args:
        displacements: Output of :func:`ltdb_levy.preprocess.compute_displacements`.
        params: Segmentation thresholds.

    Returns:
        ``(runs, audit)``. Every run appears in ``runs`` with an ``included``
        flag and an ``exclusion_reason``; nothing is dropped silently.
    """
    audit = SegmentAudit()
    angle_threshold = params.turning_angle_rad
    pieces: List[pd.DataFrame] = []

    for frag_uid, grp in displacements.groupby("fragment_uid", sort=False):
        grp = grp.sort_values("step_index", kind="mergesort")
        length = grp["length_um"].to_numpy(dtype=float)
        speed = grp["speed_um_s"].to_numpy(dtype=float)
        angles = grp["turn_angle_rad"].to_numpy(dtype=float)
        vec = grp[["dx_um", "dy_um", "dz_um"]].to_numpy(dtype=float)
        dt = grp["dt_s"].to_numpy(dtype=float)
        frame_from = grp["frame_from"].to_numpy()
        frame_to = grp["frame_to"].to_numpy()

        moving = _moving_mask(length, speed, params.speed_threshold)
        audit.non_moving_steps += int((~moving).sum())

        spans = _split_indices(moving, angles, angle_threshold)
        if not spans:
            continue

        first = grp.iloc[0]
        rows: List[Dict[str, object]] = []
        for run_id, (s, e) in enumerate(spans):
            run_vec = vec[s : e + 1].sum(axis=0)
            run_length = float(np.linalg.norm(run_vec))
            path_length = float(length[s : e + 1].sum())
            duration = float(dt[s : e + 1].sum())
            n_steps = int(e - s + 1)
            rows.append(
                {
                    "file": first["file"],
                    "cohort": first["cohort"],
                    "video": first["video"],
                    "population": first["population"],
                    "video_id": first["video_id"],
                    "track_id": first["track_id"],
                    "track_uid": first["track_uid"],
                    "fragment_uid": frag_uid,
                    "run_id": run_id,
                    "run_uid": "{}:{}".format(frag_uid, run_id),
                    "step_start": int(s),
                    "step_end": int(e),
                    "n_steps": n_steps,
                    "frame_from": int(frame_from[s]),
                    "frame_to": int(frame_to[e]),
                    "run_vector_x": float(run_vec[0]),
                    "run_vector_y": float(run_vec[1]),
                    "run_vector_z": float(run_vec[2]),
                    "run_length_um": run_length,
                    "path_length_um": path_length,
                    "straightness": (run_length / path_length) if path_length > 0 else np.nan,
                    "run_duration_s": duration,
                    "mean_speed_um_s": (path_length / duration) if duration > 0 else np.nan,
                    "is_first": run_id == 0,
                    "is_last": run_id == len(spans) - 1,
                    "included": True,
                    "exclusion_reason": "",
                }
            )
        pieces.append(pd.DataFrame(rows, columns=list(RUN_COLUMNS)))

    if not pieces:
        return pd.DataFrame(columns=list(RUN_COLUMNS)), audit

    runs = pd.concat(pieces, ignore_index=True)
    audit.n_runs_total = int(len(runs))

    # -- exclusions, applied in a fixed order so the reason is deterministic --
    reason = pd.Series([""] * len(runs), index=runs.index, dtype=object)

    if params.exclude_first_and_last_run:
        # A single-run fragment is both first and last: it is fully censored.
        boundary = runs["is_first"] | runs["is_last"]
        reason = reason.mask((reason == "") & boundary, "first_or_last_run_censored")

    too_few = runs["n_steps"] < params.minimum_run_frames
    reason = reason.mask((reason == "") & too_few, "fewer_than_minimum_run_frames")

    too_short = runs["run_length_um"] < params.minimum_run_length
    reason = reason.mask((reason == "") & too_short, "shorter_than_minimum_run_length")

    runs["exclusion_reason"] = reason
    runs["included"] = reason == ""

    audit.n_runs_included = int(runs["included"].sum())
    audit.excluded_first_last = int((reason == "first_or_last_run_censored").sum())
    audit.excluded_too_few_frames = int((reason == "fewer_than_minimum_run_frames").sum())
    audit.excluded_too_short = int((reason == "shorter_than_minimum_run_length").sum())

    logger.info(
        "Segmented %d runs (%d included) at turning angle %.1f deg, speed threshold %s",
        audit.n_runs_total,
        audit.n_runs_included,
        params.turning_angle_degrees,
        params.speed_threshold,
    )
    return runs, audit


def segmentation_sweep(
    displacements: pd.DataFrame,
    turning_angles_degrees: List[float],
    speed_thresholds: List[Optional[float]],
    base_params: SegmentParams,
) -> pd.DataFrame:
    """Segment under a grid of thresholds and summarise the resulting runs.

    The relocation-length distribution is a *function* of these thresholds, so a
    conclusion that is not stable across the grid is not a conclusion. This
    produces the evidence for that check. Thresholds are never selected by
    looking at the fitted exponent.

    Returns:
        One row per (angle, speed) combination with run counts and length
        summary statistics.
    """
    rows: List[Dict[str, object]] = []
    for angle in turning_angles_degrees:
        for speed in speed_thresholds:
            params = SegmentParams(
                turning_angle_degrees=angle,
                speed_threshold=speed,
                minimum_run_frames=base_params.minimum_run_frames,
                minimum_run_length=base_params.minimum_run_length,
                exclude_first_and_last_run=base_params.exclude_first_and_last_run,
            )
            runs, audit = segment_runs(displacements, params)
            inc = runs[runs["included"]] if len(runs) else runs
            lengths = inc["run_length_um"].to_numpy() if len(inc) else np.array([])
            rows.append(
                {
                    "turning_angle_degrees": angle,
                    "speed_threshold": speed if speed is not None else np.nan,
                    "n_runs_total": audit.n_runs_total,
                    "n_runs_included": audit.n_runs_included,
                    "n_tracks": int(inc["track_uid"].nunique()) if len(inc) else 0,
                    "runs_per_track": (
                        float(len(inc) / inc["track_uid"].nunique()) if len(inc) else np.nan
                    ),
                    "length_min": float(lengths.min()) if lengths.size else np.nan,
                    "length_median": float(np.median(lengths)) if lengths.size else np.nan,
                    "length_mean": float(lengths.mean()) if lengths.size else np.nan,
                    "length_p99": float(np.percentile(lengths, 99)) if lengths.size else np.nan,
                    "length_max": float(lengths.max()) if lengths.size else np.nan,
                    "frac_above_1um": (
                        float((lengths >= 1.0).mean()) if lengths.size else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)
