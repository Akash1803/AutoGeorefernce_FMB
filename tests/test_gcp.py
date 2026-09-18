import math

import numpy as np
import pytest

from autogeoref import anchors, edges, fit, gcp, paths, sheets

VILLAGE = "35_04_077"
CRS = 'PROJCRS["WGS 84 / UTM zone 44N"]'


def test_every_corner_gets_a_point_even_with_no_image_match():
    empty = edges.SegmentIndex([], step=0.5)
    got = gcp.corner_gcps(VILLAGE, "46B", 0.0, np.array([0.0, 0.0]), empty)
    assert len(got) >= 4, "a red parcel must still have corners to drag"
    assert all(g.matched is False for g in got)
    assert all(g.source == "auto" for g in got)


def test_matched_flag_is_set_where_imagery_agrees():
    v = sheets.outline(sheets.load_sheet(VILLAGE, "46B"))
    segs = []
    for i in range(len(v)):
        a, b = v[i], v[(i + 1) % len(v)]
        segs.append(edges.Segment(a, b, math.atan2(b[1] - a[1], b[0] - a[0]) % math.pi, math.dist(a, b)))
    idx = edges.SegmentIndex(segs, step=0.5)
    got = gcp.corner_gcps(VILLAGE, "46B", 0.0, np.array([0.0, 0.0]), idx)
    assert any(g.matched for g in got)


def test_write_then_read_round_trips_through_the_qgis_format(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    rows = [gcp.Gcp("46B", 1, (0.0, 0.0), (392800.0, 1413300.0), True, 0.0, "auto"),
            gcp.Gcp("46B", 2, (30.0, 0.0), (392830.0, 1413300.0), True, 0.0, "auto")]
    gcp.write(VILLAGE, "46B", rows, CRS)
    assert paths.points_path(VILLAGE, "46B").exists()
    assert paths.gcp_csv(VILLAGE, "46B").exists()
    assert gcp.team_points(VILLAGE, "46B") is None, "the tool's own file is not a team correction"


def test_refit_uses_the_teams_points_and_stays_rigid(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    src = np.array([(0.0, 0.0), (30.0, 0.0), (30.0, 20.0)])
    dst = fit.transform_points(src, 41.0, np.array([392800.0, 1413300.0]))
    anchors.write_points(paths.points_path(VILLAGE, "46B"),
                         [(d[0], d[1], s[0], s[1]) for s, d in zip(src, dst)], CRS)
    theta, t, rms, n = gcp.refit(VILLAGE, "46B")
    assert theta == pytest.approx(41.0, abs=0.01)
    assert t == pytest.approx(np.array([392800.0, 1413300.0]), abs=0.05)
    assert rms < 0.01 and n == 3


def test_refit_refuses_a_single_point(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    anchors.write_points(paths.points_path(VILLAGE, "46B"), [(392800.0, 1413300.0, 0.0, 0.0)], CRS)
    assert gcp.refit(VILLAGE, "46B") is None, "one GCP cannot fix a pose"


def test_a_team_row_is_never_overwritten_by_a_later_auto_run(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / VILLAGE).mkdir(parents=True)
    team = [gcp.Gcp("46B", 1, (0.0, 0.0), (392800.0, 1413300.0), True, 0.0, "team")]
    gcp.write(VILLAGE, "46B", team, CRS)
    auto = [gcp.Gcp("46B", 1, (0.0, 0.0), (999999.0, 999999.0), True, 0.0, "auto")]
    gcp.write(VILLAGE, "46B", auto, CRS)
    rows = list(paths.gcp_csv(VILLAGE, "46B").read_text(encoding="utf-8").splitlines())
    assert "392800" in rows[1] and "999999" not in rows[1]
    assert rows[1].endswith("team") or ",team," in rows[1]
