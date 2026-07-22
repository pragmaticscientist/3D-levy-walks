"""Parser and inspector for LTDB ground-truth trajectory files.

Format (Pizzagalli et al., Scientific Data 5:180129 (2018), Table 8):

    line 1   VideoID ; dx ; dy ; dz ; dt      calibration
    line 2   ch0 ; ch1 ; ch2 ; ch3 ; ch4      channel-visibility flags
    line 3+  TrackID ; x ; y ; z ; t          trajectory rows

Units, per the paper: ``x, y, z`` are already in micrometres ("position of the
cell with respect to the top-up-left most corner of the z-stack"), so they must
NOT be multiplied by dx/dy/dz -- doing so would double-scale them. ``dx``/``dy``
are the voxel size, ``dz`` the z-plane spacing, and ``dt`` the frame interval in
seconds. ``t`` is an integer frame index, so physical time is ``t * dt``.

Line 2 is a flag row, not a trajectory point. Its values vary between files
(``1;0;1;0;0``, ``0;1;0;0;0``, ``0;0;0;0;0``, ...), so it can only be identified
positionally -- never by its content.

Known corpus anomalies, all recorded rather than silently repaired:
  * ``LTDB004_b_GT.csv`` declares ``LTDB005_b`` on line 1. The filename is
    authoritative; the mismatch is flagged.
  * ``SQUARE_GT..csv`` (note the doubled dot) is a 6-point synthetic fixture.
  * Tracks are stored consecutively, not interleaved; per-track start frames
    vary and tracks may overlap in time.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..errors import ParseError
from ..logging_utils import get_logger

logger = get_logger(__name__)

#: Candidate delimiters, most likely first. LTDB uses ';' but we sniff rather
#: than assume, since the spec warns that formatting may vary between files.
_CANDIDATE_DELIMITERS: Tuple[str, ...] = (";", ",", "\t", "|")

_EXPECTED_FIELDS = 5

#: Columns of the tidy observation table produced by :func:`read_dataset`.
OBSERVATION_COLUMNS: Tuple[str, ...] = (
    "file",
    "cohort",
    "video",
    "population",
    "video_id",
    "track_id",
    "track_uid",
    "x_um",
    "y_um",
    "z_um",
    "frame",
    "time_s",
)

_STEM_RE = re.compile(r"^(?P<cohort>[A-Za-z]+)(?P<number>\d*)(?:_(?P<population>[ab]))?$")


@dataclass(frozen=True)
class LtdbFileMeta:
    """Per-file calibration and provenance."""

    source_path: Path
    file_name: str
    video_id: str          # e.g. "LTDB012_a" -- derived from the FILENAME
    video: str             # e.g. "LTDB012"
    population: str        # "a", "b", or "" when the file has no suffix
    cohort: str            # "LTDB", "CS", "SQUARE"
    internal_id: str       # the identifier declared on line 1
    id_mismatch: bool      # True when line 1 disagrees with the filename
    dx: float
    dy: float
    dz: float
    dt: float
    delimiter: str
    channel_flags: Tuple[float, ...]
    n_data_rows: int


def _parse_stem(file_name: str) -> Tuple[str, str, str, str]:
    """Derive (video_id, video, population, cohort) from a file name.

    Handles ``SQUARE_GT..csv`` (doubled dot) as well as the regular
    ``LTDB012_a_GT.csv`` / ``CS001_a_GT.csv`` forms.
    """
    name = file_name
    if name.lower().endswith(".csv"):
        name = name[: -len(".csv")]
    name = name.rstrip(".")
    if name.endswith("_GT"):
        name = name[: -len("_GT")]
    name = name.rstrip(".")

    match = _STEM_RE.match(name)
    if match is None:
        # Unrecognised naming: keep the whole stem as the id and leave the
        # structured fields blank rather than guessing.
        return name, name, "", ""

    cohort = match.group("cohort")
    number = match.group("number") or ""
    population = match.group("population") or ""
    video = "{}{}".format(cohort, number)
    video_id = "{}_{}".format(video, population) if population else video
    return video_id, video, population, cohort


def _sniff_delimiter(lines: Sequence[str]) -> str:
    """Pick the delimiter that yields a consistent field count across lines.

    Prefers a delimiter producing exactly :data:`_EXPECTED_FIELDS` fields on the
    data rows; falls back to whichever gives the most consistent split.
    """
    best: Optional[Tuple[int, int, str]] = None  # (exact_matches, consistency, delim)
    for delim in _CANDIDATE_DELIMITERS:
        counts = [len(line.split(delim)) for line in lines if line.strip()]
        if not counts:
            continue
        if max(counts) < 2:
            continue  # delimiter absent
        exact = sum(1 for c in counts if c == _EXPECTED_FIELDS)
        consistency = sum(1 for c in counts if c == counts[0])
        candidate = (exact, consistency, delim)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    if best is None:
        raise ParseError("<unknown>", 1, "could not determine a delimiter")
    return best[2]


def _read_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return fh.read().splitlines()


def _to_float(token: str, path: Path, line_number: int, field: str) -> float:
    token = token.strip()
    if token == "":
        raise ParseError(str(path), line_number, "empty value for {}".format(field))
    try:
        return float(token)
    except ValueError as exc:
        raise ParseError(
            str(path), line_number, "could not parse {} from {!r}".format(field, token)
        ) from exc


def read_ltdb_file(
    path: Path,
    coordinates_are_physical: bool = True,
) -> Tuple[LtdbFileMeta, pd.DataFrame]:
    """Parse one LTDB ground-truth file.

    Args:
        path: Path to a ``*_GT.csv`` file.
        coordinates_are_physical: When True (the documented LTDB convention),
            x/y/z are taken as micrometres and left unscaled. When False they
            are multiplied by dx/dy/dz respectively; this is offered only so the
            assumption can be tested, and is logged loudly.

    Returns:
        ``(meta, observations)`` where ``observations`` has the columns listed
        in :data:`OBSERVATION_COLUMNS`, sorted by ``(track_id, frame)``.

    Raises:
        ParseError: on a malformed header, an unparsable value, or a data row
            with the wrong number of fields.
    """
    path = Path(path)
    lines = _read_lines(path)
    if len(lines) < 3:
        raise ParseError(
            str(path), len(lines), "file has fewer than 3 lines; no trajectory rows present"
        )

    delimiter = _sniff_delimiter(lines[:15])

    # -- line 1: calibration -------------------------------------------------
    header = next(csv.reader([lines[0]], delimiter=delimiter))
    if len(header) < _EXPECTED_FIELDS:
        raise ParseError(
            str(path),
            1,
            "calibration line has {} fields, expected {} (VideoID;dx;dy;dz;dt)".format(
                len(header), _EXPECTED_FIELDS
            ),
        )
    internal_id = header[0].strip()
    dx = _to_float(header[1], path, 1, "dx")
    dy = _to_float(header[2], path, 1, "dy")
    dz = _to_float(header[3], path, 1, "dz")
    dt = _to_float(header[4], path, 1, "dt")

    # -- line 2: channel-visibility flags (skipped as data) -----------------
    flag_row = next(csv.reader([lines[1]], delimiter=delimiter))
    channel_flags = tuple(
        _to_float(tok, path, 2, "channel flag") for tok in flag_row[:_EXPECTED_FIELDS]
    )

    # -- lines 3+: trajectory rows ------------------------------------------
    track_ids: List[int] = []
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    frames: List[int] = []

    for offset, raw_line in enumerate(lines[2:]):
        line_number = offset + 3
        if not raw_line.strip():
            continue  # tolerate stray blank lines; counted via n_data_rows
        fields = next(csv.reader([raw_line], delimiter=delimiter))
        if len(fields) != _EXPECTED_FIELDS:
            raise ParseError(
                str(path),
                line_number,
                "data row has {} fields, expected {} (TrackID;x;y;z;t)".format(
                    len(fields), _EXPECTED_FIELDS
                ),
            )
        track_ids.append(int(_to_float(fields[0], path, line_number, "TrackID")))
        xs.append(_to_float(fields[1], path, line_number, "x"))
        ys.append(_to_float(fields[2], path, line_number, "y"))
        zs.append(_to_float(fields[3], path, line_number, "z"))
        frames.append(int(_to_float(fields[4], path, line_number, "t")))

    video_id, video, population, cohort = _parse_stem(path.name)
    id_mismatch = internal_id != "" and internal_id != video_id

    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    z = np.asarray(zs, dtype=float)

    if not coordinates_are_physical:
        logger.warning(
            "%s: coordinates_are_physical=False -- multiplying x,y,z by dx,dy,dz. "
            "The LTDB paper states coordinates are ALREADY in micrometres; this "
            "will double-scale them.",
            path.name,
        )
        x = x * dx
        y = y * dy
        z = z * dz

    frame_arr = np.asarray(frames, dtype=np.int64)
    track_arr = np.asarray(track_ids, dtype=np.int64)

    observations = pd.DataFrame(
        {
            "file": path.name,
            "cohort": cohort,
            "video": video,
            "population": population,
            "video_id": video_id,
            "track_id": track_arr,
            "track_uid": ["{}#{}".format(video_id, t) for t in track_arr],
            "x_um": x,
            "y_um": y,
            "z_um": z,
            "frame": frame_arr,
            "time_s": frame_arr.astype(float) * dt,
        },
        columns=list(OBSERVATION_COLUMNS),
    )
    observations = observations.sort_values(["track_id", "frame"], kind="mergesort").reset_index(
        drop=True
    )

    meta = LtdbFileMeta(
        source_path=path.resolve(),
        file_name=path.name,
        video_id=video_id,
        video=video,
        population=population,
        cohort=cohort,
        internal_id=internal_id,
        id_mismatch=id_mismatch,
        dx=dx,
        dy=dy,
        dz=dz,
        dt=dt,
        delimiter=delimiter,
        channel_flags=channel_flags,
        n_data_rows=len(observations),
    )

    if id_mismatch:
        logger.warning(
            "%s: line 1 declares id %r but the filename implies %r; using the filename.",
            path.name,
            internal_id,
            video_id,
        )

    return meta, observations


def discover_files(
    directory: Path,
    pattern: str = "*_GT.csv",
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
) -> List[Path]:
    """List LTDB files in ``directory``, applying include/exclude filters.

    ``include`` and ``exclude`` match on the bare file name. ``include`` wins
    when non-empty: only listed files are kept.
    """
    directory = Path(directory)
    paths = sorted(directory.glob(pattern))
    include_set = set(include or [])
    exclude_set = set(exclude or [])

    kept: List[Path] = []
    for p in paths:
        if include_set and p.name not in include_set:
            continue
        if p.name in exclude_set:
            logger.info("Excluding %s (listed in selection.exclude_files)", p.name)
            continue
        kept.append(p)
    return kept


def inspect_file(path: Path, coordinates_are_physical: bool = True) -> Dict[str, object]:
    """Produce a machine-readable inspection record for one file.

    Reports structure, ranges, and data-quality signals without modifying
    anything. Parse failures are captured as a row with ``ok=False`` rather than
    raised, so one bad file cannot abort the inspection of a whole corpus.
    """
    path = Path(path)
    record: Dict[str, object] = {"file": path.name, "path": str(path)}
    try:
        meta, obs = read_ltdb_file(path, coordinates_are_physical=coordinates_are_physical)
    except (ParseError, OSError, ValueError, StopIteration) as exc:
        record.update({"ok": False, "error": str(exc)})
        return record

    # Per-track diagnostics: non-monotone time, frame gaps, duplicate frames.
    non_monotone = 0
    frame_gaps = 0
    duplicate_frames = 0
    for _, grp in obs.groupby("track_id", sort=False):
        f = grp["frame"].to_numpy()
        diffs = np.diff(f)
        non_monotone += int(np.sum(diffs < 0))
        duplicate_frames += int(np.sum(diffs == 0))
        frame_gaps += int(np.sum(diffs > 1))

    coords = obs[["x_um", "y_um", "z_um"]].to_numpy()

    record.update(
        {
            "ok": True,
            "error": "",
            "delimiter": repr(meta.delimiter),
            "internal_id": meta.internal_id,
            "video_id": meta.video_id,
            "cohort": meta.cohort,
            "video": meta.video,
            "population": meta.population,
            "id_mismatch": meta.id_mismatch,
            "dx_um": meta.dx,
            "dy_um": meta.dy,
            "dz_um": meta.dz,
            "dt_s": meta.dt,
            "channel_flags": ";".join("{:g}".format(v) for v in meta.channel_flags),
            "n_observations": int(len(obs)),
            "n_tracks": int(obs["track_id"].nunique()),
            "x_min": float(coords[:, 0].min()) if len(obs) else float("nan"),
            "x_max": float(coords[:, 0].max()) if len(obs) else float("nan"),
            "y_min": float(coords[:, 1].min()) if len(obs) else float("nan"),
            "y_max": float(coords[:, 1].max()) if len(obs) else float("nan"),
            "z_min": float(coords[:, 2].min()) if len(obs) else float("nan"),
            "z_max": float(coords[:, 2].max()) if len(obs) else float("nan"),
            "frame_min": int(obs["frame"].min()) if len(obs) else -1,
            "frame_max": int(obs["frame"].max()) if len(obs) else -1,
            "n_missing_values": int(obs[["x_um", "y_um", "z_um"]].isna().to_numpy().sum()),
            "n_duplicated_rows": int(obs.duplicated(subset=["track_id", "frame"]).sum()),
            "n_duplicate_frames_within_track": duplicate_frames,
            "n_non_monotone_times": non_monotone,
            "n_frame_gaps": frame_gaps,
            # Evidence for the units convention: if z were a voxel index rather
            # than micrometres, successive z values would differ by ~1 instead
            # of by ~dz. Reported so the assumption stays auditable.
            "z_step_median": _median_abs_step(obs, "z_um"),
            "z_steps_look_like_dz": _steps_match_spacing(obs, "z_um", meta.dz),
        }
    )
    return record


def _median_abs_step(obs: pd.DataFrame, column: str) -> float:
    steps: List[np.ndarray] = []
    for _, grp in obs.groupby("track_id", sort=False):
        v = grp[column].to_numpy()
        if v.size > 1:
            d = np.abs(np.diff(v))
            steps.append(d[d > 0])
    if not steps:
        return float("nan")
    allsteps = np.concatenate(steps)
    return float(np.median(allsteps)) if allsteps.size else float("nan")


def _steps_match_spacing(obs: pd.DataFrame, column: str, spacing: float, tol: float = 1e-6) -> bool:
    """True when non-zero steps in ``column`` are integer multiples of ``spacing``.

    This is the positive evidence that the coordinate is already expressed in
    micrometres (slice position) rather than as a slice index.
    """
    if spacing <= 0:
        return False
    for _, grp in obs.groupby("track_id", sort=False):
        v = grp[column].to_numpy()
        if v.size < 2:
            continue
        d = np.abs(np.diff(v))
        d = d[d > 0]
        if d.size == 0:
            continue
        ratio = d / spacing
        if np.all(np.abs(ratio - np.round(ratio)) < tol):
            return True
    return False


def read_dataset(
    directory: Path,
    pattern: str = "*_GT.csv",
    include: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    coordinates_are_physical: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Read every selected LTDB file into one tidy observation table.

    Returns:
        ``(observations, file_table)``. ``file_table`` carries one row per file
        with its calibration and provenance, including the ``id_mismatch`` flag.
    """
    paths = discover_files(directory, pattern, include, exclude)
    if not paths:
        raise FileNotFoundError(
            "No files matching {!r} in {} after include/exclude filters".format(pattern, directory)
        )

    frames: List[pd.DataFrame] = []
    metas: List[Dict[str, object]] = []
    for p in paths:
        meta, obs = read_ltdb_file(p, coordinates_are_physical=coordinates_are_physical)
        frames.append(obs)
        metas.append(
            {
                "file": meta.file_name,
                "path": str(meta.source_path),
                "cohort": meta.cohort,
                "video": meta.video,
                "population": meta.population,
                "video_id": meta.video_id,
                "internal_id": meta.internal_id,
                "id_mismatch": meta.id_mismatch,
                "dx_um": meta.dx,
                "dy_um": meta.dy,
                "dz_um": meta.dz,
                "dt_s": meta.dt,
                "delimiter": meta.delimiter,
                "n_observations": meta.n_data_rows,
                "n_tracks": int(obs["track_id"].nunique()),
            }
        )

    observations = pd.concat(frames, ignore_index=True)
    file_table = pd.DataFrame(metas)
    logger.info(
        "Read %d files, %d observations, %d tracks",
        len(paths),
        len(observations),
        observations["track_uid"].nunique(),
    )
    return observations, file_table
