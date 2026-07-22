"""Core analysis figures generated only from versioned pipeline artifacts."""

from __future__ import annotations

import os
import tempfile
import warnings
from pathlib import Path
from typing import List

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ltdb_levy_matplotlib")
)
warnings.filterwarnings(
    "ignore", message="Unable to import Axes3D", category=UserWarning
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import Config
from .dataio import ArtifactStore


def _finish(fig: plt.Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def _condition_label(value: str) -> str:
    """Compact a verbose auditable condition key for plot labels."""
    fields = {}
    for component in str(value).split(" | "):
        if "=" in component:
            key, field_value = component.split("=", 1)
            fields[key] = field_value
    if not fields:
        return str(value)
    organ = {
        "popliteal lymph node": "PLN",
        "spleen": "spleen",
    }.get(fields.get("organ", ""), fields.get("organ", ""))
    stimulus = {
        "HIV-infected humanized T cell": "HIV T-cell",
        "Influenza Vaccine": "influenza",
        "Ovalbumin": "ovalbumin",
        "Steady State": "steady",
        "Vaccinia Virus": "vaccinia",
    }.get(fields.get("stimulus", ""), fields.get("stimulus", ""))
    dt = fields.get("dt_seconds", "")
    return "{} · {} · {} · {}s".format(
        fields.get("cell_type", ""),
        organ,
        stimulus,
        dt,
    )


def generate_figures(config: Config, store: ArtifactStore) -> List[Path]:
    """Generate the currently supported fit, pool, and pilot figures."""
    directory = store.stage_dir("report") / "figures"
    outputs: List[Path] = []

    config_hash = config.hash()
    fits = store.read_frame(
        "fit",
        "mixture_fits",
        schema_name="fits",
        expected_config_hash=config_hash,
    )
    fits = fits[fits["status"] == "ok"].sort_values("mu_hat")
    if not fits.empty:
        fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.32 * len(fits))))
        y = np.arange(len(fits))
        low = fits["mu_hat"] - fits["mu_ci_low"]
        high = fits["mu_ci_high"] - fits["mu_hat"]
        ax.errorbar(
            fits["mu_hat"],
            y,
            xerr=np.vstack((low, high)),
            fmt="o",
            capsize=2,
            color="#2b6cb0",
        )
        ax.axvline(2.0, color="#c53030", linestyle="--", linewidth=1, label="μ = 2")
        ax.set_yticks(y)
        ax.set_yticklabels([_condition_label(value) for value in fits["condition"]])
        ax.set_xlabel("Effective fitted exponent μ")
        ax.set_title("C-mixture effective exponents (track-clustered 95% CIs)")
        ax.legend(frameon=False)
        outputs.append(_finish(fig, directory / "effective_mu_estimates.png"))

        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ordered = fits.sort_values("mixture_minus_best_alternative_log_score")
        colors = [
            "#c53030" if value < 0 else "#2f855a"
            for value in ordered["mixture_minus_best_alternative_log_score"]
        ]
        ax.bar(
            np.arange(len(ordered)),
            ordered["mixture_minus_best_alternative_log_score"],
            color=colors,
        )
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_xticks(np.arange(len(ordered)))
        ax.set_xticklabels(
            [_condition_label(value) for value in ordered["condition"]],
            rotation=40,
            ha="right",
        )
        ax.set_ylabel("Mixture − best alternative\nheld-out mean log density")
        ax.set_title("Track-held-out predictive comparison")
        outputs.append(_finish(fig, directory / "model_comparison.png"))

    sensitivity = store.read_frame(
        "fit", "crossover_sensitivity", expected_config_hash=config_hash
    )
    if not sensitivity.empty:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for condition, group in sensitivity[sensitivity["status"] == "ok"].groupby(
            "condition", sort=True
        ):
            group = group.sort_values("crossover")
            ax.plot(
                group["crossover"],
                group["mu_hat"],
                marker="o",
                markersize=2.5,
                linewidth=0.9,
                alpha=0.75,
                label=_condition_label(condition),
            )
        ax.set_xscale("log", base=2)
        ax.set_xlabel("Imposed crossover (µm)")
        ax.set_ylabel("Fitted effective μ")
        ax.set_title("Sensitivity to the imposed crossover")
        ax.legend(fontsize=6, ncol=2, frameon=False)
        outputs.append(_finish(fig, directory / "crossover_sensitivity.png"))

    segmentation_fit_path = store.frame_path(
        "fit", "segmentation_fit_sensitivity"
    )
    if segmentation_fit_path.is_file():
        segmentation_fit = store.read_frame(
            "fit",
            "segmentation_fit_sensitivity",
            expected_config_hash=config_hash,
        )
        segmentation_fit = segmentation_fit[
            segmentation_fit["status"] == "ok"
        ].copy()
        if not segmentation_fit.empty:
            segmentation_fit["grid_label"] = segmentation_fit.apply(
                lambda row: "{}° / {}".format(
                    "{:g}".format(row["turning_angle_degrees"]),
                    "none"
                    if pd.isna(row["speed_threshold"])
                    else "≥{:g}".format(row["speed_threshold"]),
                ),
                axis=1,
            )
            labels = list(
                dict.fromkeys(segmentation_fit["grid_label"].tolist())
            )
            x_map = {label: index for index, label in enumerate(labels)}
            fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True)
            for condition, group in segmentation_fit.groupby(
                "condition", sort=True
            ):
                group = group.copy()
                group["x"] = group["grid_label"].map(x_map)
                group = group.sort_values("x")
                label = _condition_label(condition)
                axes[0].plot(
                    group["x"],
                    group["mu_hat"],
                    marker="o",
                    linestyle="none",
                    label=label,
                )
                axes[1].plot(
                    group["x"],
                    group["ks_bootstrap_p"],
                    marker="o",
                    linestyle="none",
                    label=label,
                )
            axes[0].axhline(
                2.0, color="#c53030", linestyle="--", linewidth=1
            )
            axes[0].set_ylabel("Effective μ")
            axes[1].axhline(
                0.05, color="#c53030", linestyle="--", linewidth=1
            )
            axes[1].set_yscale("log")
            axes[1].set_ylabel("Refitted KS bootstrap p")
            axes[1].set_xlabel("Turning angle / speed cutoff (µm s⁻¹)")
            axes[1].set_xticks(np.arange(len(labels)))
            axes[1].set_xticklabels(labels, rotation=35, ha="right")
            axes[0].legend(frameon=False, fontsize=6, ncol=2)
            fig.suptitle("Segmentation sensitivity of mixture diagnostics")
            outputs.append(
                _finish(
                    fig,
                    directory / "segmentation_fit_sensitivity.png",
                )
            )

    segmentation = store.read_frame(
        "preprocess",
        "segmentation_sweep",
        expected_config_hash=config_hash,
    )
    if not segmentation.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
        speed_values = sorted(
            segmentation["speed_threshold"].dropna().unique().tolist()
        )
        for speed in [None] + speed_values:
            subset = (
                segmentation[segmentation["speed_threshold"].isna()]
                if speed is None
                else segmentation[segmentation["speed_threshold"] == speed]
            ).sort_values("turning_angle_degrees")
            label = "no speed cutoff" if speed is None else "speed ≥ {:g}".format(speed)
            axes[0].plot(
                subset["turning_angle_degrees"],
                subset["n_runs_included"],
                marker="o",
                label=label,
            )
            axes[1].plot(
                subset["turning_angle_degrees"],
                subset["length_median"],
                marker="o",
                label=label,
            )
        axes[0].set_ylabel("Included directed runs")
        axes[1].set_ylabel("Median run length (µm)")
        for ax in axes:
            ax.set_xlabel("Turning-angle threshold (degrees)")
        axes[1].legend(frameon=False, fontsize=8)
        fig.suptitle("Directed-run definition sensitivity")
        outputs.append(
            _finish(fig, directory / "segmentation_sensitivity.png")
        )

    pools = store.read_frame(
        "preprocess", "pools", expected_config_hash=config_hash
    )
    if not pools.empty:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))
        for condition, group in pools[pools["pool_eligible"].astype(bool)].groupby(
            "condition", sort=True
        ):
            values = np.sort(group["run_length_um"].to_numpy(dtype=float))
            survival = (len(values) - np.arange(len(values))) / len(values)
            ax.step(values, survival, where="post", linewidth=0.8, alpha=0.7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Directed-run end-to-end length (µm)")
        ax.set_ylabel("Empirical CCDF")
        ax.set_title("Eligible condition relocation-length pools")
        outputs.append(_finish(fig, directory / "pool_ccdf.png"))

    runs = store.read_frame(
        "preprocess",
        "runs",
        schema_name="runs",
        expected_config_hash=config_hash,
    )
    included = runs[runs["included"].astype(bool)]
    vectors = included[["run_vector_x", "run_vector_y", "run_vector_z"]].to_numpy(
        dtype=float
    )
    norms = np.linalg.norm(vectors, axis=1)
    cosine_z = vectors[norms > 0.0, 2] / norms[norms > 0.0]
    if cosine_z.size:
        fig, ax = plt.subplots(figsize=(7.0, 4.6))
        ax.hist(cosine_z, bins=30, density=True, alpha=0.75, color="#2b6cb0")
        ax.axhline(0.5, color="#c53030", linestyle="--", label="isotropic expectation")
        ax.set_xlabel("cosine with z axis")
        ax.set_ylabel("Density")
        ax.set_title("Empirical directed-run slab anisotropy")
        ax.legend(frameon=False)
        outputs.append(_finish(fig, directory / "direction_slab_anisotropy.png"))

    replay_stage = (
        "replay"
        if store.frame_path("replay", "global_outcome_summaries").is_file()
        else "replay_pilot"
    )
    replay = store.read_frame(
        replay_stage,
        "global_outcome_summaries",
        expected_config_hash=config_hash,
    )
    if not replay.empty:
        shapes = list(config.replay.shapes)
        fig, axes = plt.subplots(
            1,
            len(shapes),
            figsize=(5.2 * len(shapes), 4.6),
            sharey=True,
        )
        if len(shapes) == 1:
            axes = [axes]
        for ax, shape in zip(axes, shapes):
            subset = replay[replay["shape"] == shape]
            for trajectory_class, group in subset.groupby(
                "trajectory_class", sort=True
            ):
                group = group.sort_values("projected_area")
                ax.plot(
                    group["projected_area"],
                    group["detection_probability"],
                    marker="o",
                    linewidth=1,
                    markersize=3,
                    label=trajectory_class,
                )
            ax.set_xscale("log", base=2)
            ax.set_title(shape)
            ax.set_xlabel("Face-on projected area")
            ax.set_ylim(bottom=-0.002)
        axes[0].set_ylabel("Detection probability")
        axes[-1].legend(fontsize=7, frameon=False, loc="best")
        fig.suptitle(
            "{} replay detection rates".format(
                "Confirmatory" if replay_stage == "replay" else "Pilot"
            )
        )
        outputs.append(_finish(fig, directory / "replay_detection_probability.png"))

    return outputs
