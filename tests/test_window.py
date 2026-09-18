import pytest

from autogeoref import anchors, fit, sheets, window

VILLAGE = "35_04_077"


def test_puvi_lookup_prefers_the_exact_letter_form_key():
    geom, kind = window.puvi_polygon(VILLAGE, "47B")
    assert geom is not None and kind == "exact"
    assert geom.area > 1000


def test_missing_puvi_polygon_is_reported_not_raised():
    geom, kind = window.puvi_polygon(VILLAGE, "99999")
    assert geom is None and kind == ""


def test_start_pose_lands_inside_the_puvi_window():
    poses, win = window.start_poses(VILLAGE, "48A")
    assert 1 <= len(poses) <= 2 and win is not None
    theta, t = poses[0]
    placed = sheets.dissolve([(None, fit.apply_pose(g, theta, t))
                              for _p, g in sheets.load_sheet(VILLAGE, "48A")])
    assert win.contains(placed.centroid)


def test_start_pose_heading_is_close_for_a_distinctive_parcel():
    a = anchors.load_anchors(VILLAGE, tool_written=set())["48A"]
    poses, _win = window.start_poses(VILLAGE, "48A")
    diffs = [abs(((theta - a.theta) + 180) % 360 - 180) for theta, _t in poses]
    assert min(diffs) < 10.0, "48A is not symmetric; the Puvi shape should find its heading"


def test_symmetric_railway_strip_keeps_both_twins():
    poses, _win = window.start_poses(VILLAGE, "170")
    assert len(poses) == 2, "a near-symmetric strip must keep the 180 degree twin"
    d = abs(((poses[0][0] - poses[1][0]) + 180) % 360 - 180)
    assert d == pytest.approx(180.0, abs=1.5)


def test_puvi_distance_is_reported_for_a_placed_parcel():
    a = anchors.load_anchors(VILLAGE, tool_written=set())["170"]
    placed = anchors.placed_geometry(VILLAGE, a)
    d = window.puvi_distance(VILLAGE, "170", placed)
    assert d is not None and 20 < d < 60, "Puvi 170 sits about 37 m from the truth"


def test_puvi_is_dropped_when_it_describes_a_different_parcel():
    """Thirukatchur 569B: the sheet is 0.98 acre, the Puvi polygon is the 259 acre village tank."""
    ok, ratio = window.puvi_trusted("35_04_074", "569B", sheet_area=3986.0)
    assert ok is False and ratio > 100
    poses, win = window.start_poses("35_04_074", "569B")
    assert poses == [] and win is None, "an untrusted Puvi polygon must not seed the search"


def test_puvi_is_kept_when_the_areas_agree():
    ok, ratio = window.puvi_trusted("35_04_074", "526A", sheet_area=81730.0)
    assert ok is True and 0.9 < ratio < 1.1
