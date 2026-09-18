import math

import numpy as np
from shapely.geometry import Polygon

from autogeoref import sheets


def test_outline_is_clockwise_and_keeps_every_vertex():
    # a square with an extra vertex in the middle of the north edge (nearly collinear)
    poly = Polygon([(0, 0), (10, 0), (10, 10), (5, 10.02), (0, 10)])
    v = sheets.outline([(None, poly)])
    assert len(v) == 5, "no collinear merge: the 10.02 vertex must survive"
    area = 0.5 * sum(v[i][0] * v[(i + 1) % 5][1] - v[(i + 1) % 5][0] * v[i][1] for i in range(5))
    assert area < 0, "clockwise rings have negative shoelace area"


def test_edge_lengths_and_turn_angles_on_a_square():
    v = sheets.outline([(None, Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]))])
    assert [round(x, 6) for x in sheets.edge_lengths(v)] == [10.0] * 4
    assert [round(x, 3) for x in sheets.turn_angles(v)] == [90.0] * 4


def test_sample_outline_spacing_and_bearings():
    v = sheets.outline([(None, Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]))])
    pts, brg = sheets.sample_outline(v, step=1.0)
    assert len(pts) == 40 and len(brg) == 40
    assert np.allclose(np.linalg.norm(np.diff(pts[:10], axis=0), axis=1), 1.0)
    assert len({round(float(b) % math.pi, 6) for b in brg}) == 2, "a square has two edge directions"


def test_load_sheet_returns_properties_and_polygons():
    got = sheets.load_sheet("35_04_077", "47B")
    assert len(got) == 3
    props, geom = got[0]
    assert props["survey_no"] == "47B"
    assert geom.is_valid and geom.area > 0


def test_dissolve_merges_touching_polygons():
    a = Polygon([(0, 0), (5, 0), (5, 5), (0, 5)])
    b = Polygon([(5, 0), (10, 0), (10, 5), (5, 5)])
    assert sheets.dissolve([(None, a), (None, b)]).area == 50
