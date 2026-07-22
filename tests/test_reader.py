from __future__ import annotations

from ltdb_levy.dataio.ltdb_reader import read_ltdb_file

from conftest import REPOSITORY_ROOT


def test_square_fixture_skips_channel_row_and_parses_semicolons():
    path = REPOSITORY_ROOT / "datasets" / "GT_TRACKS" / "SQUARE_GT..csv"
    meta, observations = read_ltdb_file(path)
    assert meta.delimiter == ";"
    assert meta.video_id == "SQUARE"
    assert len(observations) == 6
    assert observations["track_uid"].nunique() == 1


def test_known_filename_header_mismatch_is_audited():
    path = REPOSITORY_ROOT / "datasets" / "GT_TRACKS" / "LTDB004_b_GT.csv"
    meta, observations = read_ltdb_file(path)
    assert meta.video_id == "LTDB004_b"
    assert meta.internal_id == "LTDB005_b"
    assert meta.id_mismatch
    assert set(observations["video_id"]) == {"LTDB004_b"}

