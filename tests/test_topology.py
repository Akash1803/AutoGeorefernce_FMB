import pytest
from shapely.geometry import Polygon
from shapely.ops import unary_union

from autogeoref import topology

A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
B_OVERLAP = Polygon([(9.7, 0), (20, 0), (20, 10), (9.7, 10)])      # 3 sqm overlap with A
B_GAP = Polygon([(10.3, 0), (20, 0), (20, 10), (10.3, 10)])        # 3 sqm gap from A


def test_check_counts_overlaps_and_gaps():
    got = topology.check({"A": A, "B": B_OVERLAP})
    assert got["overlap_pairs"] == 1
    assert got["overlap_sqm"] == pytest.approx(3.0, abs=0.01)


def test_fix_closes_a_sub_metre_overlap_between_two_free_parcels():
    out, rep = topology.fix({"A": A, "B": B_OVERLAP}, movable={"A", "B"}, rail=set())
    assert topology.check(out)["overlap_sqm"] < 0.02
    assert rep["clips"]


def test_fix_closes_a_sub_metre_gap():
    out, _rep = topology.fix({"A": A, "B": B_GAP}, movable={"A", "B"}, rail=set())
    assert out["A"].distance(out["B"]) < 0.01


def test_an_anchor_is_never_modified():
    out, _rep = topology.fix({"A": A, "B": B_OVERLAP}, movable={"B"}, rail=set())
    assert out["A"].equals(A), "A is not movable, so it must come back untouched"


def test_railway_land_keeps_its_shape_and_the_tool_placed_parcel_yields():
    out, rep = topology.fix({"169": A, "42B": B_OVERLAP}, movable={"169", "42B"}, rail={"169"})
    assert out["169"].area == pytest.approx(A.area, abs=0.01), "railway land keeps its shape"
    assert rep["rail_conflicts"], "the conflict is still reported for the team"
    assert rep["rail_conflicts"][0][2] == pytest.approx(3.0, abs=0.5)
    assert out["42B"].intersection(out["169"]).area < 0.02, "the ordinary parcel is clipped to the strip"


def test_a_team_parcel_next_to_railway_land_is_only_reported():
    out, rep = topology.fix({"169": A, "42B": B_OVERLAP}, movable={"169"}, rail={"169"})
    assert out["42B"].equals(B_OVERLAP) and out["169"].equals(A)
    assert rep["rail_conflicts"]


def test_fix_is_idempotent():
    once, _r1 = topology.fix({"A": A, "B": B_OVERLAP}, movable={"A", "B"}, rail=set())
    twice, _r2 = topology.fix(once, movable={"A", "B"}, rail=set())
    assert topology.check(twice)["overlap_sqm"] < 0.02
    assert abs(twice["B"].area - once["B"].area) < 0.01


def test_a_tiny_plot_is_never_swallowed():
    """The 2026-09-18 16:20 run snapped at 0.30 m and destroyed two sub-metre plots of 526A."""
    sliver = Polygon([(10.0, 4.0), (10.6, 4.0), (10.6, 4.7), (10.0, 4.7)])   # 0.42 sqm, like 526A/174
    out, rep = topology.fix({"A": A, "B": B_OVERLAP, "S": sliver},
                            movable={"A", "B", "S"}, rail=set())
    assert not out["S"].is_empty, "a 0.42 sqm plot must survive the fix"
    assert out["S"].area > 0.30


def test_an_edit_that_would_eat_a_polygon_is_refused_and_reported():
    small = Polygon([(9.0, 0.0), (10.0, 0.0), (10.0, 10.0), (9.0, 10.0)])    # sits inside A
    out, rep = topology.fix({"A": A, "S": small}, movable={"S"}, rail=set())
    assert not out["S"].is_empty
    assert rep["refused"], "taking most of a polygon must be refused, not done quietly"


def test_apply_to_parts_carries_a_clip_to_the_plot_that_held_the_overlap():
    a1 = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])          # west half of A
    a2 = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])        # east half, overlaps B by 0.3 m
    original = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    fixed = Polygon([(0, 0), (9.7, 0), (9.7, 10), (0, 10)])   # after clipping against B
    out = topology.apply_to_parts([({"poly_id": 1}, a1), ({"poly_id": 2}, a2)], fixed, original)
    assert out[0][2] == "" and out[0][1].equals(a1), "the west plot was not involved"
    assert out[1][2] == "clipped" and out[1][1].area == pytest.approx(47.0, abs=0.01)


def test_apply_to_parts_gives_a_filled_gap_to_the_bordering_plot():
    a1 = Polygon([(0, 0), (5, 0), (5, 10), (0, 10)])
    a2 = Polygon([(5, 0), (10, 0), (10, 10), (5, 10)])
    original = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    fixed = Polygon([(0, 0), (10.3, 0), (10.3, 10), (0, 10)])  # a 0.3 m strip filled on the east
    out = topology.apply_to_parts([({"poly_id": 1}, a1), ({"poly_id": 2}, a2)], fixed, original)
    assert out[1][2] == "filled" and out[1][1].area == pytest.approx(53.0, abs=0.01)
    assert out[0][1].equals(a1)


def test_clean_merges_a_hairline_separated_strip_and_leaves_no_arc():
    """The dual line and the curve of 2026-09-20: a filled strip 1 mm from its plot, round-capped."""
    plot = Polygon([(0, 0), (30, 0), (30, 12), (0, 12)])
    strip = Polygon([(30.001, 0), (30.4, 0), (30.4, 12), (30.001, 12)]).buffer(0.0)      # hairline gap
    round_cap = Polygon([(30.0, 12.0), (30.4, 12.0), (30.4, 12.4)]).buffer(0.3)        # an arc on top
    dirty = unary_union([plot, strip, round_cap])
    got = topology.clean(unary_union([plot, strip]))
    assert got.geom_type == "Polygon", "one ring, not a plot plus a sliver"
    assert len(got.exterior.coords) <= 6, "a rectangle plus a strip is still four corners"
    assert got.area == pytest.approx(30.4 * 12, abs=0.2)
    got2 = topology.clean(dirty)
    assert got2.geom_type == "Polygon"
    assert len(got2.exterior.coords) < len(dirty.exterior.coords), "arc vertices are thinned"
    assert got2.area == pytest.approx(dirty.area, rel=0.01)


def test_clean_keeps_sharp_corners_and_untouched_shape():
    tri = Polygon([(0, 0), (40, 0), (5, 25)])
    got = topology.clean(tri)
    assert got.area == pytest.approx(tri.area, rel=0.005)
    assert len(got.exterior.coords) == 4


def test_gap_fill_strips_have_straight_ends():
    from shapely.geometry import Polygon as P
    a = P([(0, 0), (10, 0), (10, 10), (0, 10)]); b = P([(10.3, 0), (20, 0), (20, 10), (10.3, 10)])
    out, rep = topology.fix({"A": a, "B": b}, movable={"A", "B"}, rail=set())
    filled = out["A"] if out["A"].area > a.area + 0.1 else out["B"]
    assert len(filled.exterior.coords) == 5, "a rectangle plus a straight strip is still a rectangle"
    assert filled.area == pytest.approx(100.0 + 3.0, abs=0.05)
