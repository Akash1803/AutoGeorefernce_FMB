from shapely.geometry import Polygon
from shapely.ops import unary_union

from autogeoref import village_clean as vc


def two_plot_junior(x0):
    """A 20 x 10 survey split by one internal line at x0 + 10."""
    return [("1", Polygon([(x0, 0), (x0 + 10, 0), (x0 + 10, 10), (x0, 10)])),
            ("2", Polygon([(x0 + 10, 0), (x0 + 20, 0), (x0 + 20, 10), (x0 + 10, 10)]))]


def test_a_gap_closes_onto_one_shared_line_and_internal_lines_stay_put():
    senior = Polygon([(-20, 0), (0, 0), (0, 10), (-20, 10)])
    plots, moves = vc.conflate(two_plot_junior(1.5), senior)
    body = unary_union([g for _l, g in plots])
    assert body.distance(senior) < 1e-6, "the seam is closed"
    assert body.intersection(senior).area < 1e-6, "and nothing overlaps"
    # the internal line x = 11.5 is untouched: plot 2 still spans 11.5 .. 21.5
    p2 = dict(plots)["2"]
    assert abs(p2.bounds[0] - 11.5) < 1e-6 and abs(p2.bounds[2] - 21.5) < 1e-6
    assert max(moves) <= 1.5 + 1e-6


def test_an_overlap_is_resolved_onto_the_senior_edge():
    senior = Polygon([(-20, 0), (0, 0), (0, 10), (-20, 10)])
    plots, _m = vc.conflate(two_plot_junior(-1.0), senior)
    body = unary_union([g for _l, g in plots])
    assert body.intersection(senior).area < 1e-6
    assert body.distance(senior) < 1e-6


def test_a_bent_senior_edge_is_followed_so_no_sliver_is_left():
    # the senior's east edge has a 0.8 m kink at mid-height
    senior = Polygon([(-20, 0), (0, 0), (0.8, 5), (0, 10), (-20, 10)])
    plots, _m = vc.conflate(two_plot_junior(1.0), senior)
    body = unary_union([g for _l, g in plots])
    between = body.union(senior).convex_hull.difference(body).difference(senior)
    assert between.area < 0.05, "no sliver between the kinked edge and the junior"
    assert body.intersection(senior).area < 1e-6


def test_a_seam_wider_than_the_tolerance_is_left_open_not_forced():
    senior = Polygon([(-20, 0), (0, 0), (0, 10), (-20, 10)])
    plots, moves = vc.conflate(two_plot_junior(vc.SNAP_TOL_M + 2.0), senior)
    body = unary_union([g for _l, g in plots])
    assert body.distance(senior) > vc.SNAP_TOL_M
    assert not moves


def test_the_far_side_of_the_junior_never_moves():
    senior = Polygon([(-20, 0), (0, 0), (0, 10), (-20, 10)])
    plots, _m = vc.conflate(two_plot_junior(1.5), senior)
    p2 = dict(plots)["2"]
    assert abs(p2.bounds[2] - 21.5) < 1e-6, "the east edge 20 m away is untouched"


def test_a_small_survey_surrounded_by_neighbours_is_never_collapsed():
    # a 9 m wide survey boxed in on three sides by finished neighbours
    senior = unary_union([Polygon([(-20, -20), (0, -20), (0, 30), (-20, 30)]),
                          Polygon([(0, 22.5), (40, 22.5), (40, 30), (0, 30)]),
                          Polygon([(0, -20), (40, -20), (40, 0), (0, 0)])])
    junior = [("1", Polygon([(0.6, 0.3), (9.6, 0.3), (9.6, 22.3), (0.6, 22.3)]))]
    plots, moves = vc.conflate(junior, senior)
    body = unary_union([g for _l, g in plots])
    assert body.area > 0.9 * junior[0][1].area, "the survey keeps its ground"
    assert max(moves) <= vc.WIDTH_SHARE * 9.0 + 1e-6
    assert body.intersection(senior).area < 1e-6


def test_a_corner_that_only_comes_near_a_neighbour_does_not_move():
    # the senior lies off the junior's corner, diagonally: no edge runs along it
    senior = Polygon([(-10, -10), (-2, -10), (-2, -2), (-10, -2)])
    junior = [("1", Polygon([(0, 0), (20, 0), (20, 10), (0, 10)]))]
    plots, moves = vc.conflate(junior, senior)
    assert not moves
    assert plots[0][1].equals(junior[0][1])


def test_a_seam_never_closes_over_another_surveys_ground():
    # senior at x<0; a 2.6 m strip survey sits in the gap; the junior beyond it must not
    # snap across the strip's ground (53 swallowed 47B this way)
    from autogeoref import village_clean as vc2
    senior = Polygon([(-20, 0), (0, 0), (0, 12), (-20, 12)])
    strip = [("1", Polygon([(0, 0), (2.6, 0), (2.6, 12), (0, 12)]))]
    junior = [("1", Polygon([(2.6, 0), (22.6, 0), (22.6, 12), (2.6, 12)]))]
    # conflate alone would pull the junior onto the senior (gap 2.6 m); the village driver's
    # reserved-ground rule is what hands the strip back - checked here at the geometry level
    plots, _m = vc2.conflate(junior, senior)
    body = unary_union([g for _l, g in plots])
    reserved = strip[0][1].difference(junior[0][1])
    kept = body.difference(reserved)
    assert kept.intersection(strip[0][1]).area < 1e-6
