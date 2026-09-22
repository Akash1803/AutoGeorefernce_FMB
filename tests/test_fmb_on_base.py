import numpy as np
from shapely import affinity
from shapely.geometry import Polygon

from autogeoref import fmb_on_base


def _L(scale=1.0):
    """An L shape: it fits its twin only one way round, so the turn is unambiguous."""
    return Polygon([(0, 0), (60 * scale, 0), (60 * scale, 20 * scale),
                    (20 * scale, 20 * scale), (20 * scale, 50 * scale), (0, 50 * scale)])


def test_a_sheet_is_placed_on_its_base_parcel_without_changing_size():
    sheet = _L()
    base = affinity.translate(affinity.rotate(sheet, 35.0, origin="centroid"), 400.0, -250.0)
    pose = fmb_on_base.place(sheet, base)
    placed = fmb_on_base.apply_pose(sheet, pose, base)
    assert pose["match"] > 0.97 and abs(pose["area_ratio"] - 1.0) < 0.01
    assert abs(placed.area - sheet.area) < 0.5, "printed lengths survive: no scaling"
    assert placed.centroid.distance(base.centroid) < 0.5


def test_a_sheet_describing_different_ground_is_recognised_not_forced():
    sheet = _L()
    base = affinity.scale(_L(), xfact=2.0, yfact=2.0, origin="centroid")   # four times the area
    pose = fmb_on_base.place(sheet, base)
    assert pose["area_ratio"] < fmb_on_base.AREA_BAND[0] or pose["match"] < fmb_on_base.MATCH_MIN


def test_a_parcel_that_fits_two_ways_is_flagged_rather_than_guessed():
    # a rectangle looks the same turned 180 degrees
    rect = Polygon([(0, 0), (80, 0), (80, 20), (0, 20)])
    base = affinity.translate(rect, 100.0, 100.0)
    assert fmb_on_base.place(rect, base)["ambiguous"] is True
    # the L does not
    assert fmb_on_base.place(_L(), affinity.translate(_L(), 100.0, 100.0))["ambiguous"] is False


def test_the_first_corner_of_a_parcel_is_not_lost_to_the_closing_point():
    square = Polygon([(0, 0), (50, 0), (50, 50), (0, 50)])
    found = [tuple(np.round(p, 6)) for p in fmb_on_base.corners(square)]
    assert (0.0, 0.0) in found and len(found) == 4


def test_corners_are_the_vertices_a_person_could_point_at():
    # a square has four corners; a vertex in the middle of a straight edge is not one
    square = Polygon([(0, 0), (50, 0), (50, 25), (50, 50), (0, 50)])
    c = fmb_on_base.corners(square)
    assert len(c) == 4
    assert not any(abs(p[0] - 50) < 1e-6 and abs(p[1] - 25) < 1e-6 for p in c)


def test_control_points_pair_each_base_corner_with_one_sheet_corner():
    sheet = _L()
    base = affinity.translate(sheet, 2.0, 1.0)           # a couple of metres out
    pairs = fmb_on_base.pair_corners(sheet, base)
    assert len(pairs) == len(fmb_on_base.corners(base))
    assert all(d <= fmb_on_base.GCP_REACH_M for _s, _b, d in pairs)
    seen = [tuple(np.round(b, 6)) for _s, b, _d in pairs]
    assert len(set(seen)) == len(seen), "a base corner is used once"


def test_corners_further_apart_than_the_reach_are_not_paired():
    sheet = _L()
    base = affinity.translate(sheet, 40.0, 33.0)     # not a multiple of the shape's own steps
    assert fmb_on_base.pair_corners(sheet, base) == []


def test_the_thresholds_are_the_ones_the_sop_states():
    assert fmb_on_base.MATCH_MIN == 0.60
    assert fmb_on_base.AREA_BAND == (0.80, 1.25)
    assert fmb_on_base.GCP_REACH_M == 5.0
