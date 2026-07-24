"""Command-line interface for the fixed-arrest-scale experiment."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Dict, Sequence

from .arrest_scale_experiment import (
    load_arrest_scale_settings,
    make_arrest_scale_figures,
    run_arrest_scale_experiment,
)
from .errors import LtdbLevyError


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ltdb-arrest-scale",
        description="Run or plot the selected fixed-arrest-scale replay",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser(
        "run",
        help="Run the complete experiment from preprocessed LTDB artifacts",
    )
    run.add_argument(
        "--config",
        type=Path,
        default=Path("configs/arrest_scale.yaml"),
        help="Arrest-scale YAML settings file",
    )
    run.add_argument("--output-root", type=Path, default=None)
    run.add_argument("--finite-replays", type=int, default=None)
    run.add_argument("--batch-size", type=int, default=None)
    run.add_argument("--workers", type=int, default=None)
    run.add_argument("--clustered-bootstraps", type=int, default=None)
    run.add_argument("--gof-bootstraps", type=int, default=None)
    run.add_argument("--seed", type=int, default=None)
    run.add_argument(
        "--no-resume",
        action="store_true",
        help="Ignore a matching completion marker and rerun the experiment",
    )

    plot = subparsers.add_parser(
        "plot",
        help="Regenerate figures from completed experiment tables",
    )
    plot.add_argument(
        "--config",
        type=Path,
        default=Path("configs/arrest_scale.yaml"),
        help="Arrest-scale YAML settings file",
    )
    plot.add_argument("--output-root", type=Path, default=None)
    return parser


def _run_overrides(args: argparse.Namespace) -> Dict[str, object]:
    names = {
        "output_root": "output_root",
        "finite_replays": "finite_replays_per_track",
        "batch_size": "finite_batch_size",
        "workers": "workers",
        "clustered_bootstraps": "clustered_bootstraps",
        "gof_bootstraps": "goodness_of_fit_bootstraps",
        "seed": "random_seed",
    }
    overrides = {
        setting_name: getattr(args, argument_name)
        for argument_name, setting_name in names.items()
        if getattr(args, argument_name) is not None
    }
    if args.no_resume:
        overrides["resume"] = False
    return overrides


def main(argv: Sequence[str] = None) -> int:
    """Run the CLI and return a process exit status."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        settings = load_arrest_scale_settings(args.config)
        if args.command == "run":
            settings = replace(settings, **_run_overrides(args))
            settings.validate()
            result = run_arrest_scale_experiment(settings)
        else:
            output_root = args.output_root or settings.output_root
            paths = make_arrest_scale_figures(output_root)
            result = {
                "status": "complete",
                "output_root": str(output_root.resolve()),
                "figures": [str(path) for path in paths],
            }
    except (LtdbLevyError, OSError, ValueError) as exc:
        parser.exit(2, "error: {}\n".format(exc))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
