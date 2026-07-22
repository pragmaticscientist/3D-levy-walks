from __future__ import annotations

import numpy as np
import pandas as pd

from ltdb_levy.c_compare import ingest_c_outputs, simulate_unbounded_mixture
from ltdb_levy.targets import Target


def test_unbounded_mixture_returns_finite_first_passage_results():
    rng = np.random.default_rng(912)
    target = Target("Ball", dimension=8.0, side_length=10.0)

    result = simulate_unbounded_mixture(
        target,
        mu=2.0,
        lmax=4.0,
        n_trials=200,
        rng=rng,
        chunk_steps=8,
        max_chunks=200,
    )

    assert result.detection_distance.shape == (200,)
    assert np.all(np.isfinite(result.detection_distance))
    assert np.all(result.detection_distance >= 0.0)
    assert np.all(result.first_hit_step >= 1)


def test_c_output_ingestion_summarises_trials_and_audits_missing(tmp_path):
    present = tmp_path / "present.csv"
    pd.DataFrame(
        {
            "TargetShape": ["Ball", "Ball", "Disk", "Disk"],
            "surface": [64.0, 64.0, 64.0, 64.0],
            "detection_time": [10.0, 14.0, 20.0, 24.0],
        }
    ).to_csv(present, index=False)
    index = pd.DataFrame(
        [
            {
                "condition": "present",
                "mu_hat": 1.8,
                "c_lmax": 20,
                "num_runs": 2,
                "output_path": present.name,
                "analysis_hash": "abc",
            },
            {
                "condition": "missing",
                "mu_hat": 2.0,
                "c_lmax": 10,
                "num_runs": 2,
                "output_path": "missing.csv",
                "analysis_hash": "abc",
            },
        ]
    )
    index_path = tmp_path / "index.csv"
    index.to_csv(index_path, index=False)

    summary, audit = ingest_c_outputs(
        index_path, tmp_path, expected_analysis_hash="abc"
    )

    assert len(summary) == 2
    ball = summary[summary["shape"] == "Ball"].iloc[0]
    assert ball["mean_detection_distance"] == 12.0
    assert ball["n_trials"] == 2
    assert dict(zip(audit["condition"], audit["status"])) == {
        "present": "ok",
        "missing": "missing",
    }
