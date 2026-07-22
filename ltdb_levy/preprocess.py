"""Track cleaning and frame-to-frame displacement computation.

Responsibilities, in order:

1. Sort each track by time and resolve duplicated timestamps.
2. Split a track into *fragments* wherever consecutive observations are more
   than one frame apart. A missing-frame gap must never be interpreted as one
   long relocation -- that would manufacture spurious heavy-tail events, which
   is precisely the artefact this study must avoid.
3. Compute per-step displacement vectors, lengths, elapsed times, speeds, and
   turning angles between consecutive non-zero displacements.

Every observation or track that is dropped is counted into an audit table; the
pipeline never silently discards data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd

from .errors import ConfigError
from .logging_utils import get_logger

logger = get_logger(__name__)

DISPLACEMENT_COLUMNS: Tuple[str, ...] = (
    "file",
    "cohort",
    "video",
    "population",
    "video_id",
    "track_id",
    "track_uid",
    "fragment_id",
    "fragment_uid",
    "step_index",
    "frame_from",
    "frame_to",
    "time_s",
    "dt_s",
    "dx_um",
    "dy_um",
    "dz_um",
    "length_um",
    "speed_um_s",
    "turn_angle_rad",
)

_ID_COLUMNS = ("file", "cohort", "video", "population", "video_id", "track_id", "track_uid")


@dataclass
class PreprocessAudit:
    """Counts of every exclusion made during preprocessing."""

    duplicate_rows_dropped: int = 0
    duplicate_frames_resolved: int = 0
    tracks_below_min_observations: int = 0
    observations_in_dropped_tracks: int = 0
    fragments_created_by_gaps: int = 0
    zero_length_steps: int = 0
    rows: List[dict] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.rows is None:
            self.rows = []

    def record(self, **kwargs: object) -> None:
        self.rows.append(dict(kwargs))

    def to_frame(self) -> pd.DataFrame:
        summary = pd.DataFrame(
            [
                {"metric": "duplicate_rows_dropped", "count": self.duplicate_rows_dropped},
                {"metric": "duplicate_frames_resolved", "count": self.duplicate_frames_resolved},
                {
                    "metric": "tracks_below_min_observations",
                    "count": self.tracks_below_min_observations,
                },
                {
                    "metric": "observations_in_dropped_tracks",
                    "count": self.observations_in_dropped_tracks,
                },
                {"metric": "fragments_created_by_gaps", "count": self.fragments_created_by_gaps},
                {"metric": "zero_length_steps", "count": self.zero_length_steps},
            ]
        )
        return summary


def resolve_duplicates(
    observations: pd.DataFrame,
    policy: str = "error",
    audit: PreprocessAudit = None,
) -> pd.DataFrame:
    """Resolve duplicated ``(track_uid, frame)`` rows.

    Args:
        observations: Tidy observation table.
        policy: ``"error"`` (default) raises; ``"drop"`` keeps the first
            occurrence; ``"mean"`` averages the coordinates.
        audit: Optional audit accumulator.

    Raises:
        ConfigError: when ``policy == "error"`` and duplicates are present.
    """
    audit = audit if audit is not None else PreprocessAudit()
    dup_mask = observations.duplicated(subset=["track_uid", "frame"], keep=False)
    n_dup = int(dup_mask.sum())
    if n_dup == 0:
        return observations

    offending = (
        observations.loc[dup_mask, ["file", "track_uid", "frame"]]
        .drop_duplicates()
        .head(10)
        .to_dict("records")
    )
    if policy == "error":
        raise ConfigError(
            "Found {} duplicated (track, frame) observations. Set "
            "preprocessing.duplicate_policy to 'drop' or 'mean' to proceed. "
            "First offenders: {}".format(n_dup, offending)
        )

    if policy == "drop":
        out = observations.drop_duplicates(subset=["track_uid", "frame"], keep="first")
        audit.duplicate_rows_dropped += int(len(observations) - len(out))
    elif policy == "mean":
        group_keys = ["track_uid", "frame"]
        agg = {
            c: "first"
            for c in observations.columns
            if c not in group_keys and c not in ("x_um", "y_um", "z_um")
        }
        agg.update({"x_um": "mean", "y_um": "mean", "z_um": "mean"})
        out = (
            observations.groupby(group_keys, as_index=False, sort=False)
            .agg(agg)
            .reindex(columns=observations.columns)
        )
        audit.duplicate_rows_dropped += int(len(observations) - len(out))
    else:  # pragma: no cover - validated upstream
        raise ConfigError("Unknown duplicate_policy: {!r}".format(policy))

    audit.duplicate_frames_resolved += n_dup
    logger.warning("Resolved %d duplicated (track, frame) rows with policy %r", n_dup, policy)
    return out.reset_index(drop=True)


def split_at_gaps(
    observations: pd.DataFrame,
    split_at_missing_frames: bool = True,
    max_gap_frames: int = 1,
    audit: PreprocessAudit = None,
) -> pd.DataFrame:
    """Assign a ``fragment_id`` that increments at every frame gap.

    A gap of more than ``max_gap_frames`` between consecutive observations of a
    track starts a new fragment. With ``split_at_missing_frames=False`` the
    whole track becomes a single fragment (offered only for sensitivity
    analysis -- it bridges gaps and inflates long relocations).
    """
    audit = audit if audit is not None else PreprocessAudit()
    obs = observations.sort_values(["track_uid", "frame"], kind="mergesort").reset_index(drop=True)
    if obs.empty:
        obs["fragment_id"] = pd.Series(dtype=np.int64)
        obs["fragment_uid"] = pd.Series(dtype=object)
        return obs

    if not split_at_missing_frames:
        obs["fragment_id"] = 0
    else:
        frame = obs["frame"].to_numpy()
        track = obs["track_uid"].to_numpy()
        new_track = np.empty(len(obs), dtype=bool)
        new_track[0] = True
        new_track[1:] = track[1:] != track[:-1]

        gap = np.zeros(len(obs), dtype=bool)
        gap[1:] = (frame[1:] - frame[:-1]) > max_gap_frames
        gap &= ~new_track  # a track boundary is not a gap

        # fragment_id restarts at 0 for each track and increments on each gap
        frag = np.zeros(len(obs), dtype=np.int64)
        counter = 0
        for i in range(len(obs)):
            if new_track[i]:
                counter = 0
            elif gap[i]:
                counter += 1
            frag[i] = counter
        obs["fragment_id"] = frag
        audit.fragments_created_by_gaps += int(gap.sum())

    obs["fragment_uid"] = obs["track_uid"].astype(str) + "/" + obs["fragment_id"].astype(str)
    return obs


def filter_short_fragments(
    observations: pd.DataFrame,
    minimum_observations: int,
    audit: PreprocessAudit = None,
) -> pd.DataFrame:
    """Drop fragments with fewer than ``minimum_observations`` points.

    Exclusions are counted and recorded per fragment in the audit.
    """
    audit = audit if audit is not None else PreprocessAudit()
    sizes = observations.groupby("fragment_uid")["frame"].transform("size")
    keep = sizes >= minimum_observations
    dropped = observations.loc[~keep]

    if len(dropped):
        for frag_uid, grp in dropped.groupby("fragment_uid", sort=False):
            audit.record(
                stage="preprocess",
                fragment_uid=frag_uid,
                track_uid=grp["track_uid"].iloc[0],
                file=grp["file"].iloc[0],
                n_observations=int(len(grp)),
                reason="fewer than {} observations".format(minimum_observations),
            )
        audit.tracks_below_min_observations += int(dropped["fragment_uid"].nunique())
        audit.observations_in_dropped_tracks += int(len(dropped))
        logger.info(
            "Dropped %d fragments (%d observations) below the %d-observation minimum",
            dropped["fragment_uid"].nunique(),
            len(dropped),
            minimum_observations,
        )
    return observations.loc[keep].reset_index(drop=True)


def turning_angles(vectors: np.ndarray) -> np.ndarray:
    """Angles between consecutive 3D vectors, in radians.

    Uses ``atan2(||u x v||, u . v)`` rather than ``arccos(u.v/|u||v|)``. The
    arccos form loses catastrophic precision near 0 and pi -- exactly where a
    segmentation threshold lives -- because its derivative diverges there.

    Args:
        vectors: ``(n, 3)`` array of displacement vectors.

    Returns:
        ``(n,)`` array where element ``i`` is the angle between ``vectors[i-1]``
        and ``vectors[i]``. Element 0 is NaN. Any angle involving a zero-length
        vector is NaN, since it is undefined.
    """
    n = len(vectors)
    out = np.full(n, np.nan, dtype=float)
    if n < 2:
        return out

    u = vectors[:-1]
    v = vectors[1:]
    cross = np.cross(u, v)
    cross_norm = np.linalg.norm(cross, axis=1)
    dot = np.einsum("ij,ij->i", u, v)
    angles = np.arctan2(cross_norm, dot)

    # Undefined when either vector has zero length.
    undefined = (np.linalg.norm(u, axis=1) == 0.0) | (np.linalg.norm(v, axis=1) == 0.0)
    angles[undefined] = np.nan

    out[1:] = angles
    return out


def compute_displacements(
    observations: pd.DataFrame,
    audit: PreprocessAudit = None,
) -> pd.DataFrame:
    """Compute per-step displacements within each fragment.

    Returns a table with one row per consecutive observation pair, carrying the
    displacement vector, its length, elapsed time, speed, and the turning angle
    relative to the previous step of the same fragment.
    """
    audit = audit if audit is not None else PreprocessAudit()
    pieces: List[pd.DataFrame] = []

    for frag_uid, grp in observations.groupby("fragment_uid", sort=False):
        grp = grp.sort_values("frame", kind="mergesort")
        if len(grp) < 2:
            continue

        pos = grp[["x_um", "y_um", "z_um"]].to_numpy(dtype=float)
        frames = grp["frame"].to_numpy()
        times = grp["time_s"].to_numpy(dtype=float)

        vec = np.diff(pos, axis=0)
        length = np.linalg.norm(vec, axis=1)
        dt = np.diff(times)
        with np.errstate(divide="ignore", invalid="ignore"):
            speed = np.where(dt > 0, length / dt, np.nan)

        angles = turning_angles(vec)

        first = grp.iloc[0]
        piece = pd.DataFrame(
            {
                "file": first["file"],
                "cohort": first["cohort"],
                "video": first["video"],
                "population": first["population"],
                "video_id": first["video_id"],
                "track_id": first["track_id"],
                "track_uid": first["track_uid"],
                "fragment_id": first["fragment_id"],
                "fragment_uid": frag_uid,
                "step_index": np.arange(len(vec), dtype=np.int64),
                "frame_from": frames[:-1],
                "frame_to": frames[1:],
                "time_s": times[:-1],
                "dt_s": dt,
                "dx_um": vec[:, 0],
                "dy_um": vec[:, 1],
                "dz_um": vec[:, 2],
                "length_um": length,
                "speed_um_s": speed,
                "turn_angle_rad": angles,
            },
            columns=list(DISPLACEMENT_COLUMNS),
        )
        pieces.append(piece)

    if not pieces:
        return pd.DataFrame(columns=list(DISPLACEMENT_COLUMNS))

    displacements = pd.concat(pieces, ignore_index=True)
    audit.zero_length_steps += int((displacements["length_um"] == 0.0).sum())
    logger.info(
        "Computed %d displacements across %d fragments (%d zero-length)",
        len(displacements),
        displacements["fragment_uid"].nunique(),
        int((displacements["length_um"] == 0.0).sum()),
    )
    return displacements


def preprocess(
    observations: pd.DataFrame,
    duplicate_policy: str = "error",
    split_at_missing_frames: bool = True,
    minimum_observations_per_track: int = 6,
) -> Tuple[pd.DataFrame, pd.DataFrame, PreprocessAudit]:
    """Run the full preprocessing chain.

    Returns:
        ``(cleaned_observations, displacements, audit)``.
    """
    audit = PreprocessAudit()
    obs = resolve_duplicates(observations, policy=duplicate_policy, audit=audit)
    obs = split_at_gaps(obs, split_at_missing_frames=split_at_missing_frames, audit=audit)
    obs = filter_short_fragments(obs, minimum_observations_per_track, audit=audit)
    displacements = compute_displacements(obs, audit=audit)
    return obs, displacements, audit
