from __future__ import annotations

import pandas as pd
import pytest

from ltdb_levy.dataio.artifacts import ArtifactError, ArtifactStore


def test_artifact_roundtrip_and_integrity_check(tmp_path):
    store = ArtifactStore(tmp_path)
    original = pd.DataFrame({"x": [1, 2], "label": ["a", "b"]})
    path = store.write_frame("stage", "observations", original, "observations")
    loaded = store.read_frame("stage", "observations", "observations")
    pd.testing.assert_frame_equal(loaded, original)

    with path.open("a", encoding="utf-8") as fh:
        fh.write("3,c\n")
    with pytest.raises(ArtifactError, match="hash mismatch"):
        store.read_frame("stage", "observations", "observations")


def test_incremental_artifact_writer_commits_bounded_chunks(tmp_path):
    store = ArtifactStore(tmp_path)
    writer = store.incremental_frame_writer(
        "stage", "replay_trials", schema_name="replay_trials"
    )
    writer.append(pd.DataFrame({"trial": [1, 2], "value": [3.0, 4.0]}))
    writer.append(pd.DataFrame({"trial": [3], "value": [5.0]}))
    path = writer.close()
    assert path.is_file()
    loaded = store.read_frame("stage", "replay_trials", "replay_trials")
    assert loaded["trial"].tolist() == [1, 2, 3]


def test_run_context_rejects_stale_stage(tmp_path):
    store = ArtifactStore(tmp_path)
    store.write_json("stage", "run_context", {"config_hash": "old"})
    with pytest.raises(ArtifactError, match="current hash"):
        store.require_run_context("stage", "new")
