import pytest
from shapely.geometry import Polygon

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


def test_railway_overlap_is_listed_not_clipped():
    out, rep = topology.fix({"169": A, "42B": B_OVERLAP}, movable={"169", "42B"}, rail={"169"})
    assert out["169"].area == pytest.approx(A.area, abs=0.01), "railway land keeps its shape"
    assert rep["rail_conflicts"], "the conflict is reported for the team"
    assert rep["rail_conflicts"][0][2] == pytest.approx(3.0, abs=0.5)


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
