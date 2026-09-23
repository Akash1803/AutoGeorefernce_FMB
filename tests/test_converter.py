import math

import fmb_to_geojson as F


def _seg(p0, p1, layer="subdivision_line"):
    return {"p0": p0, "p1": p1, "layer": layer, "lw": 1.0, "dashed": False}


def test_drop_arrows_removes_an_offset_arrow_and_keeps_the_plot_lines():
    tol = 2.7
    square = [_seg((0, 0), (100, 0), "survey_boundary_main"), _seg((100, 0), (100, 100), "survey_boundary_main"),
              _seg((100, 100), (0, 100), "survey_boundary_main"), _seg((0, 100), (0, 0), "survey_boundary_main")]
    divider = [_seg((50, 0), (50, 100))]
    # arrow: foot inside the square at (48, 60), tip outside at (-20, 60), two head strokes leaning
    # back along the line at about 23 degrees on opposite sides, and the bar closing the head
    arrow = [_seg((48, 60), (-20, 60)), _seg((-20, 60), (-13, 63)), _seg((-20, 60), (-13, 57)), _seg((-13, 63), (-13, 57))]
    kept, dropped = F.drop_arrows(square + divider + arrow, tol)
    assert len(dropped) == 4 and len(kept) == 5
    assert all(k in square + divider for k in kept)


def test_drop_arrows_keeps_a_line_with_a_single_tick():
    tol = 2.7
    lines = [_seg((0, 0), (0, 100)), _seg((0, 100), (5, 95))]
    kept, dropped = F.drop_arrows(lines, tol)
    assert dropped == [] and len(kept) == 2


def test_merge_slivers_folds_a_thin_face_into_its_longest_neighbour():
    big = {"area": 100.0, "rings": [[[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]]}
    other = {"area": 100.0, "rings": [[[10.5, 0], [20.5, 0], [20.5, 10], [10.5, 10], [10.5, 0]]]}
    sliver = {"area": 5.0, "rings": [[[10, 0], [10.5, 0], [10.5, 10], [10, 10], [10, 0]]]}   # 0.5 m wide strip
    out, merged = F.merge_slivers([big, other, sliver], min_area=2.5)
    assert merged == 0, "5 m2 is above the threshold"
    out, merged = F.merge_slivers([big, other, sliver], min_area=6.0)
    assert merged == 1 and len(out) == 2
    assert max(p["area"] for p in out) == 105.0
    assert all(abs(p["area"] - round(p["area"])) < 1e-6 for p in out)


def test_merge_slivers_leaves_an_isolated_sliver():
    lone = {"area": 1.0, "rings": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]]}
    far = {"area": 100.0, "rings": [[[50, 50], [60, 50], [60, 60], [50, 60], [50, 50]]]}
    out, merged = F.merge_slivers([lone, far], min_area=2.5)
    assert merged == 0 and len(out) == 2
