import math

import numpy as np
import pytest
from shapely.geometry import Polygon

from autogeoref import align, edges, sheets


def _square_index(x0=0.0, y0=0.0, side=40.0, step=0.5):
    """Four walls of a square, as if seen in imagery."""
    segs = [edges.Segment((x0, y0), (x0 + side, y0), 0.0, side),
            edges.Segment((x0 + side, y0), (x0 + side, y0 + side), math.pi / 2, side),
            edges.Segment((x0 + side, y0 + side), (x0, y0 + side), 0.0, side),
            edges.Segment((x0, y0 + side), (x0, y0), math.pi / 2, side)]
    return edges.SegmentIndex(segs, step=step)


SHEET = [(None, Polygon([(0, 0), (40, 0), (40, 40), (0, 40)]))]


def test_share_is_one_at_the_true_pose_and_low_when_far_away():
    pts, brg = sheets.sample_outline(sheets.outline(SHEET), step=1.0)
    idx = _square_index()
    assert align.share(pts, brg, idx, 0.0, np.array([0.0, 0.0])) > 0.95
    assert align.share(pts, brg, idx, 0.0, np.array([25.0, 0.0])) < 0.60


def _index_for(poly, step=0.5):
    """Treat a polygon's own edges as the imagery that should be matched."""
    c = list(poly.exterior.coords)
    segs = [edges.Segment(c[i], c[i + 1],
                          math.atan2(c[i + 1][1] - c[i][1], c[i + 1][0] - c[i][0]) % math.pi,
                          math.dist(c[i], c[i + 1])) for i in range(len(c) - 1)]
    return edges.SegmentIndex(segs, step=step)


def test_search_recovers_a_shifted_pose():
    pts, brg = sheets.sample_outline(sheets.outline(SHEET), step=1.0)
    res = align.search(pts, brg, _square_index(), [(0.0, np.array([6.0, -5.0]))],
                       shift_m=10.0, rot_deg=10.0)
    assert np.allclose(res.t, [0.0, 0.0], atol=0.5)
    assert res.share > 0.95


def test_a_shape_with_a_re_entrant_corner_is_not_ambiguous():
    """A distinctive corner pins both rotation and position, as the walls of 48A did."""
    poly = Polygon([(0, 0), (60, 0), (60, 20), (25, 20), (25, 50), (0, 50)])
    pts, brg = sheets.sample_outline(sheets.outline([(None, poly)]), step=1.0)
    res = align.search(pts, brg, _index_for(poly), [(0.0, np.array([6.0, -5.0]))],
                       shift_m=10.0, rot_deg=10.0)
    assert np.allclose(res.t, [0.0, 0.0], atol=0.5)
    assert res.ambiguous is False and res.observable is True


def test_a_bare_rectangle_slides_along_its_own_length_and_is_flagged():
    """Sliding a plain strip along itself barely changes the score; that is the v1 failure mode."""
    poly = Polygon([(0, 0), (100, 0), (100, 20), (0, 20)])
    pts, brg = sheets.sample_outline(sheets.outline([(None, poly)]), step=1.0)
    res = align.search(pts, brg, _index_for(poly), [(0.0, np.array([6.0, -5.0]))],
                       shift_m=10.0, rot_deg=10.0)
    assert res.share > 0.95, "the true pose is still found"
    assert res.ambiguous is True, "but it must not be trusted on its own"


def test_parallel_lines_are_reported_ambiguous_not_placed():
    """Rails and road markings repeat every few metres; the score then has many equal peaks."""
    segs = [edges.Segment((0.0, y), (200.0, y), 0.0, 200.0) for y in (0.0, 3.0, 6.0, 9.0, 12.0)]
    idx = edges.SegmentIndex(segs, step=0.5)
    strip = [(None, Polygon([(0, 0), (150, 0), (150, 3), (0, 3)]))]
    pts, brg = sheets.sample_outline(sheets.outline(strip), step=1.0)
    res = align.search(pts, brg, idx, [(0.0, np.array([0.0, 0.5]))], shift_m=8.0, rot_deg=4.0)
    assert res.ambiguous is True
    assert res.alternative is not None


def test_observability_needs_two_directions():
    pts, brg = sheets.sample_outline(sheets.outline(SHEET), step=1.0)
    full = align.observability(pts, brg, _square_index(), 0.0, np.array([0.0, 0.0]))
    assert full["ok"] is True and full["cross_m"] >= 20.0

    one_way = edges.SegmentIndex([edges.Segment((0.0, 0.0), (40.0, 0.0), 0.0, 40.0),
                                  edges.Segment((0.0, 40.0), (40.0, 40.0), 0.0, 40.0)], step=0.5)
    partial = align.observability(pts, brg, one_way, 0.0, np.array([0.0, 0.0]))
    assert partial["ok"] is False, "one direction leaves the along-track position free"
    assert partial["cross_m"] < 20.0


def test_a_pose_with_no_image_support_scores_zero_and_is_not_observable():
    pts, brg = sheets.sample_outline(sheets.outline(SHEET), step=1.0)
    empty = edges.SegmentIndex([], step=0.5)
    assert align.share(pts, brg, empty, 0.0, np.array([0.0, 0.0])) == 0.0
    obs = align.observability(pts, brg, empty, 0.0, np.array([0.0, 0.0]))
    assert obs["ok"] is False and obs["primary_m"] == 0.0
