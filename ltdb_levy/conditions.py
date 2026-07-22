"""Assignment of tracks to homogeneous experimental conditions.

A *condition* is the unit that receives one fitted exponent and one empirical
relocation-length pool. It must be homogeneous in cell type, organ, immune
stimulus and microscope sampling interval.

The biological annotation is distributed with this analysis after verification
against the LTDB data descriptor (Pizzagalli et al., Scientific Data 5:180129)
and deposited SQL database. The loader still preserves the defensive workflow:

* :func:`write_metadata_template` scaffolds a CSV with one row per (video,
  population) when the curated file is absent;
* rows carry ``verified`` and ``provenance`` columns;
* unverified rows are **not** pooled across files when
  ``require_verified_metadata`` is set -- they fall back to their own
  ``(video, population)`` condition.

Nothing here ever invents an annotation. A condition built on a guess would
silently pool distinct cell types, which is the failure mode the whole design is
meant to prevent.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .errors import MetadataError
from .logging_utils import get_logger

logger = get_logger(__name__)

METADATA_COLUMNS: Sequence[str] = (
    "video_id",
    "video",
    "population",
    "cohort",
    "cell_type",
    "organ",
    "stimulus",
    "mouse_strain",
    "dt_seconds",
    "verified",
    "provenance",
)

#: Where the annotation comes from, quoted in the template so the person
#: filling it in knows exactly which tables to read.
_PROVENANCE_HINT = (
    "Pizzagalli et al. 2018 Sci Data 5:180129 - Table 3 "
    "(channel/cell population) for cell_type and Table 4 "
    "(experimental settings) for organ and stimulus."
)


def write_metadata_template(file_table: pd.DataFrame, path: Path) -> Path:
    """Write a metadata scaffold with one row per (video, population).

    ``dt_seconds`` is filled in from each file's own calibration header, which
    is authoritative and in-repo. The biological columns are left blank for a
    human to complete from the paper.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, object]] = []
    for _, r in file_table.sort_values("video_id").iterrows():
        rows.append(
            {
                "video_id": r["video_id"],
                "video": r["video"],
                "population": r["population"],
                "cohort": r["cohort"],
                "cell_type": "",
                "organ": "",
                "stimulus": "",
                "mouse_strain": "",
                "dt_seconds": r["dt_s"],
                "verified": False,
                "provenance": "",
            }
        )
    template = pd.DataFrame(rows, columns=list(METADATA_COLUMNS))
    with path.open("w", encoding="utf-8") as fh:
        fh.write("# LTDB biological metadata. Fill cell_type/organ/stimulus from:\n")
        fh.write("#   {}\n".format(_PROVENANCE_HINT))
        fh.write("# Set verified=True only for rows you have checked against the paper.\n")
        fh.write("# dt_seconds is taken from each file's own header and is authoritative.\n")
        template.to_csv(fh, index=False)
    logger.info("Wrote metadata template for %d files to %s", len(template), path)
    return path


def load_metadata(path: Optional[Path]) -> Optional[pd.DataFrame]:
    """Load the biological metadata table, if one is configured.

    Returns ``None`` when no path is given. Raises when the file exists but is
    missing required columns, since a partially-formed metadata table is more
    dangerous than none at all.
    """
    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        raise MetadataError(
            "input.metadata_file is set to {} but that file does not exist. Run "
            "`ltdb-levy inspect` to generate a template.".format(path)
        )
    md = pd.read_csv(path, comment="#")
    missing = [c for c in ("video_id", "cell_type", "organ", "stimulus", "verified") if c not in md.columns]
    if missing:
        raise MetadataError(
            "Metadata file {} is missing required columns: {}".format(path, ", ".join(missing))
        )
    md["verified"] = md["verified"].astype(str).str.strip().str.lower().isin(("true", "1", "yes"))
    duplicates = md["video_id"].duplicated(keep=False)
    if duplicates.any():
        ids = sorted(md.loc[duplicates, "video_id"].astype(str).unique())
        raise MetadataError(
            "Metadata file {} contains duplicate video_id rows: {}".format(
                path, ", ".join(ids[:10])
            )
        )
    logger.info(
        "Loaded metadata for %d files (%d verified)", len(md), int(md["verified"].sum())
    )
    return md


def assign_conditions(
    frame: pd.DataFrame,
    metadata: Optional[pd.DataFrame],
    grouping_fields: Sequence[str],
    file_table: Optional[pd.DataFrame] = None,
    require_verified: bool = True,
    fallback_to_file: bool = True,
) -> pd.DataFrame:
    """Attach a ``condition`` column to a per-run or per-track table.

    Rows whose metadata is absent or unverified fall back to a condition named
    after their own ``video_id``, so they are still analysed -- just never
    pooled with anything else on an unverified assumption.

    Args:
        frame: Table containing at least ``video_id``.
        metadata: Optional metadata table from :func:`load_metadata`.
        grouping_fields: Fields that jointly define a condition, e.g.
            ``["cell_type", "organ", "stimulus", "dt_seconds"]``.
        file_table: Per-file calibration, used to supply ``dt_seconds`` when it
            is one of the grouping fields and no metadata is present.
        require_verified: Only pool rows flagged ``verified``.
        fallback_to_file: Fall back to per-``video_id`` conditions. When False,
            unverified rows are labelled ``condition = ""`` and excluded later.

    Returns:
        A copy of ``frame`` with ``condition`` and ``condition_source`` columns.
    """
    out = frame.copy()

    # Always make dt available: it comes from the file header, not the paper.
    if file_table is not None and "dt_seconds" not in out.columns:
        dt_map = dict(zip(file_table["video_id"], file_table["dt_s"]))
        out["dt_seconds"] = out["video_id"].map(dt_map)

    if metadata is None:
        if not fallback_to_file:
            raise MetadataError(
                "No metadata file configured and conditions.fallback_to_file is False; "
                "there is no defensible way to define a condition."
            )
        out["condition"] = out["video_id"].astype(str)
        out["condition_source"] = "file_fallback"
        logger.warning(
            "No biological metadata supplied: using per-(video, population) conditions. "
            "Fill in the metadata template to pool by cell type / organ / stimulus."
        )
        return out

    md = metadata.copy()
    usable = md["verified"] if require_verified else pd.Series(True, index=md.index)
    md = md.assign(_usable=usable)

    merge_cols = ["video_id", "_usable"] + [
        c for c in grouping_fields if c in md.columns and c not in out.columns
    ]
    out = out.merge(md[merge_cols], on="video_id", how="left")
    out["_usable"] = out["_usable"].fillna(False).astype(bool)

    missing_fields = [c for c in grouping_fields if c not in out.columns]
    if missing_fields:
        raise MetadataError(
            "Grouping fields {} are not present in the metadata or the data. "
            "Available columns: {}".format(missing_fields, sorted(out.columns))
        )

    # A row is poolable only when it is verified AND every grouping field is set.
    complete = out["_usable"].copy()
    for c in grouping_fields:
        col = out[c]
        complete &= col.notna() & (col.astype(str).str.strip() != "")

    def _label(row: pd.Series) -> str:
        return " | ".join("{}={}".format(c, row[c]) for c in grouping_fields)

    condition = pd.Series("", index=out.index, dtype=object)
    if complete.any():
        condition.loc[complete] = out.loc[complete].apply(_label, axis=1)

    source = pd.Series("metadata", index=out.index, dtype=object)
    if (~complete).any():
        if fallback_to_file:
            condition.loc[~complete] = out.loc[~complete, "video_id"].astype(str)
            source.loc[~complete] = "file_fallback"
        else:
            source.loc[~complete] = "unassigned"

    out["condition"] = condition
    out["condition_source"] = source
    out = out.drop(columns=["_usable"])

    n_fallback = int((source == "file_fallback").sum())
    if n_fallback:
        logger.warning(
            "%d rows lacked verified metadata and fell back to per-file conditions",
            n_fallback,
        )
    logger.info(
        "Assigned %d distinct conditions (%d from metadata)",
        out["condition"].nunique(),
        int((source == "metadata").sum()),
    )
    return out


def condition_summary(runs: pd.DataFrame, min_tracks: int, min_runs: int) -> pd.DataFrame:
    """Summarise each condition and flag eligibility for pool construction.

    Conditions below the configured minimums are reported with a reason rather
    than dropped, so the audit trail shows exactly what was excluded and why.
    """
    included = runs[runs["included"]] if "included" in runs.columns else runs
    rows: List[Dict[str, object]] = []
    for cond, grp in included.groupby("condition", sort=True):
        n_tracks = int(grp["track_uid"].nunique())
        n_runs = int(len(grp))
        lengths = grp["run_length_um"].to_numpy()
        reasons: List[str] = []
        if n_tracks < min_tracks:
            reasons.append("fewer than {} tracks".format(min_tracks))
        if n_runs < min_runs:
            reasons.append("fewer than {} runs".format(min_runs))
        rows.append(
            {
                "condition": cond,
                "condition_source": grp["condition_source"].iloc[0]
                if "condition_source" in grp.columns
                else "",
                "n_tracks": n_tracks,
                "n_runs": n_runs,
                "n_videos": int(grp["video_id"].nunique()),
                "runs_per_track": n_runs / n_tracks if n_tracks else np.nan,
                "length_min": float(lengths.min()) if lengths.size else np.nan,
                "length_median": float(np.median(lengths)) if lengths.size else np.nan,
                "length_max": float(lengths.max()) if lengths.size else np.nan,
                "eligible": not reasons,
                "ineligible_reason": "; ".join(reasons),
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "condition",
                "condition_source",
                "n_tracks",
                "n_runs",
                "n_videos",
                "runs_per_track",
                "length_min",
                "length_median",
                "length_max",
                "eligible",
                "ineligible_reason",
            ]
        )
    summary = pd.DataFrame(rows).sort_values(
        ["eligible", "n_runs"], ascending=[False, False]
    )
    n_elig = int(summary["eligible"].sum()) if len(summary) else 0
    logger.info(
        "%d of %d conditions meet the pool minimums (>=%d tracks, >=%d runs)",
        n_elig,
        len(summary),
        min_tracks,
        min_runs,
    )
    return summary.reset_index(drop=True)
