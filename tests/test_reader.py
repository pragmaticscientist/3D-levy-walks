from __future__ import annotations

from ltdb_levy.dataio.ltdb_reader import read_ltdb_file


def test_square_fixture_skips_channel_row_and_parses_semicolons(tmp_path):
    path = tmp_path / "SQUARE_GT..csv"
    path.write_text(
        "\n".join(
            [
                "SQUARE;0.5;0.5;2;20",
                "1;0;0;0;0",
                "1;2.5;2.5;5;1",
                "1;2.5;7.5;5;2",
                "1;2.5;12.5;5;3",
                "1;2.5;17.5;5;4",
                "1;2.5;22.5;5;5",
                "1;2.5;27.5;5;6",
            ]
        ),
        encoding="utf-8",
    )
    meta, observations = read_ltdb_file(path)
    assert meta.delimiter == ";"
    assert meta.video_id == "SQUARE"
    assert len(observations) == 6
    assert observations["track_uid"].nunique() == 1


def test_known_filename_header_mismatch_is_audited(tmp_path):
    path = tmp_path / "LTDB004_b_GT.csv"
    path.write_text(
        "\n".join(
            [
                "LTDB005_b;0.8;0.8;3;15",
                "0;1;0;0;0",
                "1;372.61;430.32;38.37;61",
            ]
        ),
        encoding="utf-8",
    )
    meta, observations = read_ltdb_file(path)
    assert meta.video_id == "LTDB004_b"
    assert meta.internal_id == "LTDB005_b"
    assert meta.id_mismatch
    assert set(observations["video_id"]) == {"LTDB004_b"}
