"""YAML configuration loading, validation, and hashing.

Every methodological choice in the pipeline is expressed here rather than
hard-coded, so that the analysis is fully described by its config file. The
config is hashed and written into the run manifest.

Python 3.9 note: we use ``typing.List``/``Optional`` rather than PEP 604 unions
because dataclass field types are introspected at runtime.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import yaml

from .errors import ConfigError

SCHEMA_VERSION = 1

# Duplicate-timestamp handling. "error" is the default because a duplicated
# timestamp in a hand-curated ground-truth file usually signals a real problem.
_DUPLICATE_POLICIES = ("error", "drop", "mean")
_LMAX_METHODS = ("empirical_max", "user", "quantile")
_WEIGHTINGS = ("cell_balanced", "step_weighted")
_OVERSHOOT_POLICIES = ("discard", "truncate")
_BUDGET_TYPES = ("path_length", "fixed_common")
_SHAPES = ("Ball", "Disk", "Line")


@dataclass(frozen=True)
class InputConfig:
    tracks_directory: Path
    glob: str = "*_GT.csv"
    metadata_file: Optional[Path] = None
    # The LTDB paper (Sci Data 5:180129, Table 8) states x, y, z are already in
    # micrometres. Scaling them by dx/dy/dz would double-count. Kept
    # configurable so the assumption is visible rather than buried in code.
    coordinates_are_physical: bool = True


@dataclass(frozen=True)
class SelectionConfig:
    include_files: List[str] = field(default_factory=list)
    exclude_files: List[str] = field(default_factory=list)
    minimum_observations_per_track: int = 6


@dataclass(frozen=True)
class PreprocessingConfig:
    split_at_missing_frames: bool = True
    duplicate_policy: str = "error"
    interpolation: bool = False


@dataclass(frozen=True)
class SegmentationConfig:
    method: str = "angle_speed"
    turning_angle_degrees: float = 60.0
    speed_threshold: Optional[float] = None
    speed_threshold_units: str = "um_per_second"
    alternative_turning_angles: List[float] = field(
        default_factory=lambda: [45.0, 60.0, 90.0]
    )
    alternative_speed_thresholds: List[Optional[float]] = field(
        default_factory=lambda: [None]
    )
    exclude_first_and_last_run: bool = True
    minimum_run_frames: int = 1
    minimum_run_length: float = 0.0


@dataclass(frozen=True)
class ConditionsConfig:
    """How homogeneous conditions are formed.

    Defaults to the biological grouping from the LTDB paper's metadata tables.
    Files without verified metadata fall back to their own (video, population)
    condition rather than being pooled on a guess.
    """

    grouping_fields: List[str] = field(
        default_factory=lambda: ["cell_type", "organ", "stimulus", "dt_seconds"]
    )
    fallback_to_file: bool = True
    require_verified_metadata: bool = True


@dataclass(frozen=True)
class EmpiricalPoolConfig:
    primary_weighting: str = "cell_balanced"
    alternative_weightings: List[str] = field(default_factory=lambda: ["step_weighted"])
    use_complete_runs_only: bool = True
    minimum_tracks_per_pool: int = 15
    minimum_runs_per_pool: int = 100
    optional_block_lengths: List[int] = field(default_factory=lambda: [3, 5])
    sample_length_duration_jointly: bool = False


@dataclass(frozen=True)
class FittingConfig:
    """Parameters of the C-mixture fit.

    The fitted model is f(l) = a(mu, lmax) * min(1, (l/crossover)^-mu), matching
    the sampler in func.c. ``crossover`` supplies the fixed value or optimizer
    start. When ``estimate_crossover`` is true, one crossover is fitted per
    condition. ``lmax`` comes from the data.
    """

    lmax_method: str = "empirical_max"
    user_lmax: Optional[float] = None
    lmax_quantile: float = 0.99
    crossover: float = 1.0
    estimate_crossover: bool = False
    crossover_bounds: Optional[List[float]] = None
    crossover_sensitivity: List[float] = field(
        default_factory=lambda: [0.25, 0.5, 1.0, 2.0, 4.0]
    )
    mu_bounds: List[float] = field(default_factory=lambda: [1.001, 10.0])
    goodness_of_fit_bootstraps: int = 500
    segmentation_goodness_of_fit_bootstraps: int = 100
    clustered_bootstraps: int = 1000
    cv_folds: int = 10
    random_seed: int = 12345


@dataclass(frozen=True)
class ReplayConfig:
    """Replay experiment parameters.

    Conventions deliberately mirror the existing C engine so that a run of the
    unmodified simulator at the fitted mu is directly comparable: detection
    radius 1, cylinder disk membership, line capsule along y, face-on
    (maximum-projection) area matching.
    """

    domain_side_length: float = 512.0
    detection_radius: float = 1.0
    shapes: List[str] = field(default_factory=lambda: ["Ball", "Disk", "Line"])
    projected_areas: List[float] = field(
        default_factory=lambda: [64.0, 256.0, 1024.0, 4096.0, 16384.0]
    )
    replays_per_track_target: int = 200
    budget_type: str = "path_length"
    fixed_budget_quantile: float = 0.5
    overshoot_policy: str = "discard"
    check_initial_position: bool = False
    block_lengths: List[int] = field(default_factory=lambda: [3, 5])
    pilot_tracks_per_condition: int = 2
    pilot_replays_per_track: int = 10
    pilot_min_detected: int = 20
    pilot_min_undetected: int = 20
    random_seed: int = 12345


@dataclass(frozen=True)
class OutputConfig:
    root: Path = Path("outputs")
    figures: bool = True


@dataclass(frozen=True)
class Config:
    schema_version: int
    run_id: str
    input: InputConfig
    selection: SelectionConfig
    preprocessing: PreprocessingConfig
    segmentation: SegmentationConfig
    conditions: ConditionsConfig
    empirical_pool: EmpiricalPoolConfig
    fitting: FittingConfig
    replay: ReplayConfig
    output: OutputConfig
    source_path: Optional[Path] = None

    # -- derived output locations -------------------------------------------

    @property
    def run_dir(self) -> Path:
        return self.output.root

    def stage_dir(self, stage: str) -> Path:
        """Return (and create) the output directory for a pipeline stage."""
        d = self.output.root / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    # -- provenance ----------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable dict of the configuration."""

        def _convert(obj: Any) -> Any:
            if isinstance(obj, Path):
                return str(obj)
            if isinstance(obj, dict):
                return {k: _convert(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [_convert(v) for v in obj]
            return obj

        d = asdict(self)
        d.pop("source_path", None)
        return _convert(d)

    def hash(self) -> str:
        """Return a stable SHA-256 hash of configuration and source inputs.

        The biological metadata is part of the analysis specification, and the
        trajectory files are the source data. Including their content digests
        makes stale-artifact checks react when either changes even though the
        YAML paths themselves remain unchanged.
        """
        selected_files = []
        include = set(self.selection.include_files)
        exclude = set(self.selection.exclude_files)
        if self.input.tracks_directory.is_dir():
            for path in sorted(self.input.tracks_directory.glob(self.input.glob)):
                if include and path.name not in include:
                    continue
                if path.name in exclude or not path.is_file():
                    continue
                selected_files.append(
                    {"name": path.name, "sha256": _sha256_file(path)}
                )

        metadata_digest = None
        if (
            self.input.metadata_file is not None
            and self.input.metadata_file.is_file()
        ):
            metadata_digest = _sha256_file(self.input.metadata_file)

        payload = {
            "config": self.to_dict(),
            "source_inputs": {
                "metadata_sha256": metadata_digest,
                "track_files": selected_files,
            },
        }
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# Loading and validation
# ---------------------------------------------------------------------------


def _require_mapping(value: Any, name: str) -> Dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ConfigError("Config section '{}' must be a mapping, got {}".format(name, type(value).__name__))
    return value


def _sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Hash a configured source without importing the artifact subsystem."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _check_choice(value: str, allowed: Sequence[str], name: str) -> str:
    if value not in allowed:
        raise ConfigError(
            "{} must be one of {}, got {!r}".format(name, ", ".join(allowed), value)
        )
    return value


def _optional_path(value: Any, base: Path) -> Optional[Path]:
    if value in (None, ""):
        return None
    return _resolve(Path(value), base)


def _resolve(p: Path, base: Path) -> Path:
    return p if p.is_absolute() else (base / p)


def load_config(path: Path) -> Config:
    """Load, validate, and return the pipeline configuration.

    Relative paths inside the config are resolved against the config file's own
    directory, so a config is portable alongside its data.

    Raises:
        ConfigError: if the file is missing or any value is invalid.
    """
    path = Path(path)
    if not path.is_file():
        raise ConfigError("Configuration file not found: {}".format(path))

    try:
        with path.open("r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        raise ConfigError("Could not parse YAML in {}: {}".format(path, exc)) from exc

    if not isinstance(raw, dict):
        raise ConfigError("Configuration root must be a mapping in {}".format(path))

    base = path.parent.resolve()

    schema_version = int(raw.get("schema_version", SCHEMA_VERSION))
    if schema_version != SCHEMA_VERSION:
        raise ConfigError(
            "Unsupported schema_version {} (this build understands {})".format(
                schema_version, SCHEMA_VERSION
            )
        )

    run_id = str(raw.get("run_id", "primary"))

    # -- input --------------------------------------------------------------
    raw_input = _require_mapping(raw.get("input"), "input")
    if "tracks_directory" not in raw_input:
        raise ConfigError("input.tracks_directory is required")
    input_cfg = InputConfig(
        tracks_directory=_resolve(Path(raw_input["tracks_directory"]), base),
        glob=str(raw_input.get("glob", "*_GT.csv")),
        metadata_file=_optional_path(raw_input.get("metadata_file"), base),
        coordinates_are_physical=bool(raw_input.get("coordinates_are_physical", True)),
    )
    if not input_cfg.tracks_directory.is_dir():
        raise ConfigError(
            "input.tracks_directory does not exist: {}".format(input_cfg.tracks_directory)
        )

    # -- selection ----------------------------------------------------------
    raw_sel = _require_mapping(raw.get("selection"), "selection")
    selection_cfg = SelectionConfig(
        include_files=[str(x) for x in raw_sel.get("include_files") or []],
        exclude_files=[str(x) for x in raw_sel.get("exclude_files") or []],
        minimum_observations_per_track=int(raw_sel.get("minimum_observations_per_track", 6)),
    )
    if selection_cfg.minimum_observations_per_track < 2:
        raise ConfigError("selection.minimum_observations_per_track must be >= 2")

    # -- preprocessing ------------------------------------------------------
    raw_pre = _require_mapping(raw.get("preprocessing"), "preprocessing")
    preprocessing_cfg = PreprocessingConfig(
        split_at_missing_frames=bool(raw_pre.get("split_at_missing_frames", True)),
        duplicate_policy=_check_choice(
            str(raw_pre.get("duplicate_policy", "error")),
            _DUPLICATE_POLICIES,
            "preprocessing.duplicate_policy",
        ),
        interpolation=bool(raw_pre.get("interpolation", False)),
    )
    if preprocessing_cfg.interpolation:
        raise ConfigError(
            "preprocessing.interpolation is not implemented: interpolating across "
            "missing frames would manufacture relocations that were never observed"
        )

    # -- segmentation -------------------------------------------------------
    raw_seg = _require_mapping(raw.get("segmentation"), "segmentation")
    segmentation_cfg = SegmentationConfig(
        method=_check_choice(
            str(raw_seg.get("method", "angle_speed")), ("angle_speed",), "segmentation.method"
        ),
        turning_angle_degrees=float(raw_seg.get("turning_angle_degrees", 60.0)),
        speed_threshold=(
            None if raw_seg.get("speed_threshold") is None else float(raw_seg["speed_threshold"])
        ),
        speed_threshold_units=str(raw_seg.get("speed_threshold_units", "um_per_second")),
        alternative_turning_angles=[
            float(x) for x in raw_seg.get("alternative_turning_angles") or [45.0, 60.0, 90.0]
        ],
        alternative_speed_thresholds=[
            None if x is None else float(x)
            for x in raw_seg.get("alternative_speed_thresholds") or [None]
        ],
        exclude_first_and_last_run=bool(raw_seg.get("exclude_first_and_last_run", True)),
        minimum_run_frames=int(raw_seg.get("minimum_run_frames", 1)),
        minimum_run_length=float(raw_seg.get("minimum_run_length", 0.0)),
    )
    if not 0.0 < segmentation_cfg.turning_angle_degrees < 180.0:
        raise ConfigError("segmentation.turning_angle_degrees must be in (0, 180)")
    if any(not 0.0 < x < 180.0 for x in segmentation_cfg.alternative_turning_angles):
        raise ConfigError("segmentation.alternative_turning_angles must all be in (0, 180)")
    if segmentation_cfg.minimum_run_frames < 1:
        raise ConfigError("segmentation.minimum_run_frames must be >= 1")
    if segmentation_cfg.minimum_run_length < 0.0:
        raise ConfigError("segmentation.minimum_run_length must be >= 0")
    speed_values = [segmentation_cfg.speed_threshold] + list(
        segmentation_cfg.alternative_speed_thresholds
    )
    if any(x is not None and x < 0.0 for x in speed_values):
        raise ConfigError("segmentation speed thresholds must be non-negative")

    # -- conditions ---------------------------------------------------------
    raw_cond = _require_mapping(raw.get("conditions"), "conditions")
    conditions_cfg = ConditionsConfig(
        grouping_fields=[
            str(x)
            for x in raw_cond.get("grouping_fields")
            or ["cell_type", "organ", "stimulus", "dt_seconds"]
        ],
        fallback_to_file=bool(raw_cond.get("fallback_to_file", True)),
        require_verified_metadata=bool(raw_cond.get("require_verified_metadata", True)),
    )

    # -- empirical pool -----------------------------------------------------
    raw_pool = _require_mapping(raw.get("empirical_pool"), "empirical_pool")
    pool_cfg = EmpiricalPoolConfig(
        primary_weighting=_check_choice(
            str(raw_pool.get("primary_weighting", "cell_balanced")),
            _WEIGHTINGS,
            "empirical_pool.primary_weighting",
        ),
        alternative_weightings=[
            _check_choice(str(x), _WEIGHTINGS, "empirical_pool.alternative_weightings")
            for x in raw_pool.get("alternative_weightings") or ["step_weighted"]
        ],
        use_complete_runs_only=bool(raw_pool.get("use_complete_runs_only", True)),
        minimum_tracks_per_pool=int(raw_pool.get("minimum_tracks_per_pool", 15)),
        minimum_runs_per_pool=int(raw_pool.get("minimum_runs_per_pool", 100)),
        optional_block_lengths=[int(x) for x in raw_pool.get("optional_block_lengths") or [3, 5]],
        sample_length_duration_jointly=bool(
            raw_pool.get("sample_length_duration_jointly", False)
        ),
    )
    if pool_cfg.minimum_tracks_per_pool < 1 or pool_cfg.minimum_runs_per_pool < 1:
        raise ConfigError("empirical_pool minimum track and run counts must be positive")
    if any(x < 1 for x in pool_cfg.optional_block_lengths):
        raise ConfigError("empirical_pool.optional_block_lengths must all be >= 1")

    # -- fitting ------------------------------------------------------------
    raw_fit = _require_mapping(raw.get("fitting"), "fitting")
    fitting_cfg = FittingConfig(
        lmax_method=_check_choice(
            str(raw_fit.get("lmax_method", "empirical_max")), _LMAX_METHODS, "fitting.lmax_method"
        ),
        user_lmax=(None if raw_fit.get("user_lmax") is None else float(raw_fit["user_lmax"])),
        lmax_quantile=float(raw_fit.get("lmax_quantile", 0.99)),
        crossover=float(raw_fit.get("crossover", 1.0)),
        estimate_crossover=bool(raw_fit.get("estimate_crossover", False)),
        crossover_bounds=(
            None
            if raw_fit.get("crossover_bounds") is None
            else [float(x) for x in raw_fit["crossover_bounds"]]
        ),
        crossover_sensitivity=[
            float(x) for x in raw_fit.get("crossover_sensitivity") or [0.25, 0.5, 1.0, 2.0, 4.0]
        ],
        mu_bounds=[float(x) for x in raw_fit.get("mu_bounds") or [1.001, 10.0]],
        goodness_of_fit_bootstraps=int(raw_fit.get("goodness_of_fit_bootstraps", 500)),
        segmentation_goodness_of_fit_bootstraps=int(
            raw_fit.get("segmentation_goodness_of_fit_bootstraps", 100)
        ),
        clustered_bootstraps=int(raw_fit.get("clustered_bootstraps", 1000)),
        cv_folds=int(raw_fit.get("cv_folds", 10)),
        random_seed=int(raw_fit.get("random_seed", 12345)),
    )
    if fitting_cfg.crossover <= 0.0:
        raise ConfigError("fitting.crossover must be positive")
    if fitting_cfg.crossover_bounds is not None:
        if (
            len(fitting_cfg.crossover_bounds) != 2
            or fitting_cfg.crossover_bounds[0] <= 0.0
            or fitting_cfg.crossover_bounds[0] >= fitting_cfg.crossover_bounds[1]
        ):
            raise ConfigError(
                "fitting.crossover_bounds must be a positive increasing pair"
            )
    if len(fitting_cfg.mu_bounds) != 2 or fitting_cfg.mu_bounds[0] >= fitting_cfg.mu_bounds[1]:
        raise ConfigError("fitting.mu_bounds must be [low, high] with low < high")
    if fitting_cfg.mu_bounds[0] <= 1.0:
        # mu <= 1 makes the tail integral diverge in the mu != 1 branch; the
        # mu == 1 case is handled analytically but is not a valid bound.
        raise ConfigError("fitting.mu_bounds lower bound must exceed 1.0")
    if fitting_cfg.lmax_method == "user" and fitting_cfg.user_lmax is None:
        raise ConfigError("fitting.user_lmax is required when lmax_method is 'user'")
    if not 0.0 < fitting_cfg.lmax_quantile <= 1.0:
        raise ConfigError("fitting.lmax_quantile must be in (0, 1]")
    if (
        not fitting_cfg.estimate_crossover
        and fitting_cfg.user_lmax is not None
        and fitting_cfg.user_lmax <= fitting_cfg.crossover
    ):
        raise ConfigError("fitting.user_lmax must exceed fitting.crossover")
    if any(x <= 0.0 for x in fitting_cfg.crossover_sensitivity):
        raise ConfigError("fitting.crossover_sensitivity values must be positive")
    if (
        fitting_cfg.goodness_of_fit_bootstraps < 0
        or fitting_cfg.segmentation_goodness_of_fit_bootstraps < 0
        or fitting_cfg.clustered_bootstraps < 0
    ):
        raise ConfigError("bootstrap counts must be non-negative")

    # -- replay -------------------------------------------------------------
    raw_replay = _require_mapping(raw.get("replay"), "replay")
    replay_cfg = ReplayConfig(
        domain_side_length=float(raw_replay.get("domain_side_length", 512.0)),
        detection_radius=float(raw_replay.get("detection_radius", 1.0)),
        shapes=[
            _check_choice(str(x), _SHAPES, "replay.shapes")
            for x in raw_replay.get("shapes") or ["Ball", "Disk", "Line"]
        ],
        projected_areas=[float(x) for x in raw_replay.get("projected_areas") or []],
        replays_per_track_target=int(raw_replay.get("replays_per_track_target", 200)),
        budget_type=_check_choice(
            str(raw_replay.get("budget_type", "path_length")), _BUDGET_TYPES, "replay.budget_type"
        ),
        fixed_budget_quantile=float(raw_replay.get("fixed_budget_quantile", 0.5)),
        overshoot_policy=_check_choice(
            str(raw_replay.get("overshoot_policy", "discard")),
            _OVERSHOOT_POLICIES,
            "replay.overshoot_policy",
        ),
        check_initial_position=bool(raw_replay.get("check_initial_position", False)),
        block_lengths=[int(x) for x in raw_replay.get("block_lengths") or [3, 5]],
        pilot_tracks_per_condition=int(
            raw_replay.get("pilot_tracks_per_condition", 2)
        ),
        pilot_replays_per_track=int(raw_replay.get("pilot_replays_per_track", 10)),
        pilot_min_detected=int(raw_replay.get("pilot_min_detected", 20)),
        pilot_min_undetected=int(raw_replay.get("pilot_min_undetected", 20)),
        random_seed=int(raw_replay.get("random_seed", 12345)),
    )
    if replay_cfg.detection_radius <= 0.0:
        raise ConfigError("replay.detection_radius must be positive")
    if replay_cfg.domain_side_length <= 0.0:
        raise ConfigError("replay.domain_side_length must be positive")
    if replay_cfg.replays_per_track_target < 1:
        raise ConfigError("replay.replays_per_track_target must be >= 1")
    if replay_cfg.pilot_tracks_per_condition < 1:
        raise ConfigError("replay.pilot_tracks_per_condition must be >= 1")
    if replay_cfg.pilot_replays_per_track < 1:
        raise ConfigError("replay.pilot_replays_per_track must be >= 1")
    if replay_cfg.pilot_min_detected < 1 or replay_cfg.pilot_min_undetected < 1:
        raise ConfigError("replay pilot minimum event counts must be >= 1")

    # Every requested projected area must admit a positive raw dimension for
    # every requested shape; A > pi*d^2 is the binding constraint.
    import math

    min_area = math.pi * replay_cfg.detection_radius ** 2
    for area in replay_cfg.projected_areas:
        if area <= min_area:
            raise ConfigError(
                "replay.projected_areas contains {} which is not > pi*d^2 = {:.4f}; "
                "no target of any shape has a positive size at that area".format(area, min_area)
            )

    # -- output -------------------------------------------------------------
    raw_out = _require_mapping(raw.get("output"), "output")
    output_cfg = OutputConfig(
        root=_resolve(Path(str(raw_out.get("root", "outputs"))), base),
        figures=bool(raw_out.get("figures", True)),
    )

    return Config(
        schema_version=schema_version,
        run_id=run_id,
        input=input_cfg,
        selection=selection_cfg,
        preprocessing=preprocessing_cfg,
        segmentation=segmentation_cfg,
        conditions=conditions_cfg,
        empirical_pool=pool_cfg,
        fitting=fitting_cfg,
        replay=replay_cfg,
        output=output_cfg,
        source_path=path.resolve(),
    )
