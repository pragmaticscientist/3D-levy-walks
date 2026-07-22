from __future__ import annotations

from ltdb_levy.cli import _select_budget_quantile_tracks


def test_pilot_track_selection_spans_budget_quantiles():
    budgets = {"t{}".format(index): float(index) for index in range(8)}

    assert _select_budget_quantile_tracks(budgets, 2) == ["t2", "t6"]
    assert _select_budget_quantile_tracks(budgets, 4) == [
        "t1",
        "t3",
        "t5",
        "t7",
    ]
