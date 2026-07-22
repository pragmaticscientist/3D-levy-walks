from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from ltdb_levy.conditions import load_metadata
from ltdb_levy.config import load_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


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


def test_curated_ltdb_metadata_is_complete_and_source_verified():
    metadata = load_metadata(
        REPOSITORY_ROOT / "configs" / "ltdb_metadata.csv"
    )

    assert metadata is not None
    assert len(metadata) == 44
    assert metadata["verified"].all()
    assert not metadata[["cell_type", "organ", "stimulus"]].isna().any().any()
    by_id = metadata.set_index("video_id")
    assert by_id.loc["LTDB010_a", "cell_type"] == "B Cells"
    assert by_id.loc["LTDB012_b", "cell_type"] == "T Cells"
    assert by_id.loc["LTDB017_b", "cell_type"] == "Neutrophils"
    assert by_id.loc["CS005_a", "stimulus"] == "Vaccinia Virus"
