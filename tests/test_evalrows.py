import csv

import pytest
from shapely import affinity
from shapely.geometry import Polygon

from autogeoref import config, evalrows, references

SQ = Polygon([(0, 0), (100, 0), (100, 40), (0, 40)])


def _ref(stretched=False, disputed=False, theta=0.0):
    return references.Reference(
        village="V", survey="7A", hand_file="h.gpkg", source="points", theta=theta, tx=0.0, ty=0.0,
        rms_m=0.5, n_points=3, area_scale=1.0, stretch_long_pct=0.0, stretch_short_pct=0.0,
        extent_m=100.0, centroid_gap_m=0.0, disagree_mean_m=0.2, disagree_p95_m=0.5, disagree_max_m=0.8,
        disagree_corner_mean_m=0.3, disagree_corner_max_m=0.6, corners_matched=4, band="<=5%", stretched=stretched, stretched_reason="x" if stretched else "",
        disputed=disputed, disputed_reason="y" if disputed else "", geom_fmb=SQ, geom_hand=SQ)


STATUS = {"survey": "7A", "status": "placed", "method": "neighbour", "colour": "amber", "confidence": "60",
          "support_m": "200", "support_next_m": "50", "boundary_rms_m": "0.4", "neighbours": "1,2,3",
          "anchor_partners": "1", "side_ok": "2", "side_bad": "0", "contradictions": "0", "printed_far": "",
          "share": "0.7", "observable": "True", "ambiguous": "False", "in_window": "True",
          "sigma_pos_m": "0.3", "sigma_head_deg": "0.2", "heading_deg": "1.0", "pose_tx": "1.5", "pose_ty": "2.5",
          "n_candidates": "9", "pass": "2",
          "notes": "neighbour support 200 m, next best 50 m | topology: conformed,clipped,conform max 0.41 m"}


def test_row_carries_every_feature_and_both_errors():
    cfg = config.Config()
    placed = affinity.translate(SQ, 1.0, 0.0)
    row = evalrows.row_from_status("V", "r1", "2026-09-21T10:00:00", "seed", ["48A"], STATUS, _ref(), cfg, placed,
                                   imagery_source="google_xyz", window_sha256="abc", tool_git_sha="deadbee")
    d = row.as_dict()
    for name in ("support_m", "runner_up_ratio", "boundary_rms_m", "n_neighbours", "anchor_partners",
                 "side_ok", "side_bad", "share", "in_window", "pass_no", "err_fmb_m", "err_hand_m"):
        assert d[name] != "", name
    assert row.runner_up_ratio == pytest.approx(0.25)
    assert row.n_neighbours == 3 and row.pass_no == 2 and row.observable == 1 and row.in_window == 1
    assert row.topo_max_move_m == pytest.approx(0.41) and row.topo_clipped == 1
    assert row.err_fmb_m == pytest.approx(1.0) and row.err_hand_m == pytest.approx(1.0)
    assert row.heading_err_deg == pytest.approx(1.0)
    assert row.within_3m == 1 and row.seeds == "48A" and row.imagery_source == "google_xyz"
    assert set(evalrows.FIELDS) == set(d)


def test_label_rules():
    cfg = config.Config()
    far = affinity.translate(SQ, 5.0, 0.0)
    assert evalrows.row_from_status("V", "r", "", "loo", [], STATUS, _ref(), cfg, far).within_3m == 0
    assert evalrows.row_from_status("V", "r", "", "loo", [], STATUS, _ref(stretched=True), cfg, far).within_3m is None
    assert evalrows.row_from_status("V", "r", "", "loo", [], STATUS, _ref(disputed=True), cfg, SQ).within_3m is None
    assert evalrows.row_from_status("V", "r", "", "loo", [], STATUS, None, cfg, SQ).within_3m is None
    stretched = evalrows.row_from_status("V", "r", "", "loo", [], STATUS, _ref(stretched=True), cfg, far)
    assert stretched.err_fmb_m == pytest.approx(5.0), "errors are still reported for a stretched reference"


def test_rows_round_trip_through_the_csv(tmp_path):
    cfg = config.Config()
    row = evalrows.row_from_status("V", "r1", "", "seed", ["48A"], STATUS, _ref(), cfg, SQ)
    p = evalrows.append_rows([row], tmp_path / "rows.csv")
    evalrows.append_rows([row], p)
    back = evalrows.load_rows(p)
    assert len(back) == 2 and back[0]["survey"] == "7A" and back[0]["within_3m"] == "1"
    with p.open(encoding="utf-8") as fh:
        assert fh.readline().count("run_id") == 1, "one header only"


def test_split_is_by_village_and_folds_cover_labelled_villages_only():
    rows = [{"village": "A", "survey": "1", "within_3m": "1"}, {"village": "A", "survey": "2", "within_3m": "0"},
            {"village": "B", "survey": "1", "within_3m": "1"}, {"village": "C", "survey": "9", "within_3m": ""}]
    train, test = evalrows.split_by_village(rows, "A")
    assert {r["village"] for r in test} == {"A"} and {r["village"] for r in train} == {"B", "C"}
    fs = evalrows.folds(rows)
    assert [v for v, _tr, _te in fs] == ["A", "B"], "C has no label, so it is never a test fold"
    for v, tr, te in fs:
        assert not ({r["village"] for r in tr} & {r["village"] for r in te})


def test_backfill_reads_an_old_leave_one_out_csv(tmp_path):
    p = tmp_path / "loo.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["survey", "colour", "centroid_error_m", "heading_error_deg", "share",
                                           "observable", "boundary_rms_m", "neighbours", "method", "margin"])
        w.writeheader()
        w.writerow({"survey": "7A", "colour": "green", "centroid_error_m": "1.52", "heading_error_deg": "0.9",
                    "share": "0.71", "observable": "False", "boundary_rms_m": "0.57", "neighbours": "43A,47A,170",
                    "method": "anchors", "margin": "0.1"})
    rows = evalrows.backfill_leave_one_out(p, "V", "old", config.Config(), {"7A": _ref()})
    assert rows[0].backfilled == 1 and rows[0].err_fmb_m == pytest.approx(1.52) and rows[0].within_3m == 1
    assert rows[0].n_neighbours == 3 and rows[0].support_m is None
