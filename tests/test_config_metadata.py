from __future__ import annotations

from dataclasses import replace

from ltdb_levy.config import load_config

from conftest import REPOSITORY_ROOT


def test_analysis_hash_changes_with_metadata_and_track_content(tmp_path):
    track_directory = tmp_path / "tracks"
    track_directory.mkdir()
    track = track_directory / "TEST_a_GT.csv"
    track.write_text("source A\n", encoding="utf-8")
    metadata_path = tmp_path / "metadata.csv"
    metadata_path.write_text("metadata A\n", encoding="utf-8")

    original = load_config(REPOSITORY_ROOT / "configs" / "primary.yaml")
    config = replace(
        original,
        input=replace(
            original.input,
            tracks_directory=track_directory,
            glob="*.csv",
            metadata_file=metadata_path,
        ),
        selection=replace(
            original.selection, include_files=[], exclude_files=[]
        ),
    )
    first = config.hash()
    metadata_path.write_text("metadata B\n", encoding="utf-8")
    second = config.hash()
    track.write_text("source B\n", encoding="utf-8")
    third = config.hash()

    assert first != second
    assert second != third
