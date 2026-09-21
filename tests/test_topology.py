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
    # which side receives the strip is a tie (both border it equally), so test the result, not the side
    assert out["A"].area + out["B"].area == pytest.approx(200.0, abs=0.05), "the 3 m2 gap is filled"
    assert out["A"].distance(out["B"]) < 1e-9 and out["A"].intersection(out["B"]).area < 1e-6
    for g in out.values():
        assert len(g.exterior.coords) == 5, "a rectangle plus a straight strip is still a rectangle"



def _ring_keys(g):
    return {(round(x, 6), round(y, 6)) for x, y in g.exterior.coords}


def test_conform_makes_a_shared_boundary_identical():
    """B sits 0.3 m off A with a bend of its own; A has a mid-edge vertex B lacks."""
    A = Polygon([(0, 0), (10, 0), (10, 6), (10, 10), (0, 10)])                 # settled, vertex at (10,6)
    B = Polygon([(10.3, 0.2), (20, 0), (20, 10), (10.4, 10.3), (10.25, 4.0)])  # movable, a bend at y=4
    out, rep = topology.conform({"B": [({"poly_id": 1}, B)]}, {"A": A}, ["B"], anchors=())
    nb = out["B"][0][1]
    assert nb.is_valid and out["B"][0][2] == "conformed"
    assert nb.intersection(A).area < 1e-6 and nb.distance(A) < 1e-9, "touching, not crossing"
    keys = _ring_keys(nb)
    assert (10.0, 0.0) in keys and (10.0, 10.0) in keys, "corners snapped onto A's corners"
    assert (10.0, 6.0) in keys, "A's mid-edge vertex inserted into B"
    assert all(abs(x - 10.0) < 1e-9 for x, y in keys if x < 15), "the whole shared run lies on x = 10"
    assert rep["conformed"]["B"]["max_move_m"] < 0.6


def test_conform_keeps_a_parcels_plots_stitched():
    """Moving the outer ring must move the same vertices in every plot that shares them."""
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    left = Polygon([(10.3, 0), (15, 0), (15, 10), (10.3, 10)])
    right = Polygon([(15, 0), (20, 0), (20, 10), (15, 10)])
    out, _rep = topology.conform({"B": [({"poly_id": 1}, left), ({"poly_id": 2}, right)]}, {"A": A}, ["B"])
    nl, nr = out["B"][0][1], out["B"][1][1]
    assert nl.distance(A) < 1e-9 and nl.intersection(A).area < 1e-6
    assert nl.intersection(nr).area < 1e-6 and nl.distance(nr) < 1e-9, "internal boundary still shared"
    assert nr.equals(right), "the plot away from the shared boundary is untouched"
    from shapely.ops import unary_union
    body = unary_union([nl, nr])
    assert body.geom_type == "Polygon" and len(body.interiors) == 0, "no hole opened between the plots"


def test_conform_never_touches_the_settled_parcel_and_refuses_to_eat_a_plot():
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    sliver = Polygon([(10.3, 4.0), (10.7, 4.0), (10.7, 4.5), (10.3, 4.5)])      # 0.4 x 0.5 m plot
    big = Polygon([(10.7, 0), (20, 0), (20, 10), (10.7, 10)])
    out, rep = topology.conform({"B": [({"poly_id": 1}, sliver), ({"poly_id": 2}, big)]}, {"A": A}, ["B"])
    assert not out["B"][0][1].is_empty and out["B"][0][1].area > 0.1, "the small plot survives"
    assert A.equals(Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]))


def test_conform_orders_by_certainty_and_settles_each_parcel_for_the_next():
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    B = Polygon([(10.2, 0), (20, 0), (20, 10), (10.2, 10)])
    C = Polygon([(20.3, 0), (30, 0), (30, 10), (20.3, 10)])
    out, _rep = topology.conform({"B": [({}, B)], "C": [({}, C)]}, {"A": A}, ["B", "C"])
    nb, nc = out["B"][0][1], out["C"][0][1]
    assert nb.distance(A) < 1e-9 and nc.distance(nb) < 1e-9
    assert nb.intersection(nc).area < 1e-6 and (20.0, 0.0) in _ring_keys(nc)



def test_resolve_clips_a_real_crossing_from_the_less_certain_parcel_only():
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])                     # team parcel
    B = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])                   # certain, shares A's edge exactly
    C = Polygon([(18.5, 0), (30, 0), (30, 10), (18.5, 10)])               # 1.5 m into B: beyond tolerance
    out, rep = topology.resolve({"B": [({"poly_id": 1}, B)], "C": [({"poly_id": 1}, C)]}, {"A": A}, ["B", "C"])
    nb, nc = out["B"][0][1], out["C"][0][1]
    assert nb.equals(B), "the more certain parcel keeps its shape"
    assert nc.intersection(nb).area < 1e-6 and nc.area == pytest.approx(C.area - 15.0, abs=0.01)
    assert rep["clips"] and rep["clips"][0][0] == "C"


def test_resolve_fills_an_enclosed_sliver_by_exact_union():
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    B = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])
    # C sits below both but stops 0.3 m short between x = 4 and x = 16: a thin enclosed hole
    C = Polygon([(0, -10), (20, -10), (20, 0), (16, 0), (16, -0.3), (4, -0.3), (4, 0), (0, 0)])
    out, rep = topology.resolve({"C": [({"poly_id": 1}, C)]}, {"A": A, "B": B}, ["C"], anchors=set(), tol=0.2)
    nc = out["C"][0][1]
    from shapely.ops import unary_union
    u = unary_union([A, B, nc])
    assert u.geom_type == "Polygon" and len(u.interiors) == 0, "the sliver hole is closed"
    assert rep["fills"] and nc.area == pytest.approx(C.area + 12 * 0.3, abs=0.01)


def test_resolve_keeps_plots_of_a_parcel_exactly_stitched():
    """Nothing in the stage may make two plots of one parcel drift apart (the 46A defect)."""
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    plots = [({"poly_id": i}, Polygon([(10.3 + 2 * i, 0.1), (12.3 + 2 * i, 0.1), (12.3 + 2 * i, 10.2), (10.3 + 2 * i, 10.2)]))
             for i in range(5)]
    out, _rep = topology.resolve({"B": plots}, {"A": A}, ["B"])
    geoms = [g for _p, g, _n in out["B"]]
    from shapely.ops import unary_union
    import itertools
    assert sum(g1.intersection(g2).area for g1, g2 in itertools.combinations(geoms, 2)) < 1e-9
    body = unary_union(geoms)
    assert body.geom_type == "Polygon" and len(body.interiors) == 0
    assert body.distance(A) < 1e-9 and body.intersection(A).area < 1e-6



def test_conform_inserts_whichever_way_round_the_edge_is_keyed():
    """The mirror image of the shared-boundary case: sorted coordinate keys now run against the
    ring direction, and the edge's first vertex is one that moves. 171 lost 48A's corner this way."""
    A = Polygon([(0, 0), (-10, 0), (-10, 6), (-10, 10), (0, 10)])              # settled, west side
    B = Polygon([(-10.3, 0.2), (-20, 0), (-20, 10), (-10.4, 10.3), (-10.25, 4.0)])
    out, rep = topology.conform({"B": [({"poly_id": 1}, B)]}, {"A": A}, ["B"], anchors=())
    nb = out["B"][0][1]
    keys = _ring_keys(nb)
    assert (-10.0, 6.0) in keys, "A's mid-edge vertex inserted into B"
    assert nb.is_valid and nb.intersection(A).area < 1e-6 and nb.distance(A) < 1e-9
    assert rep["conformed"]["B"]["vertices_inserted"] >= 1



def test_conform_keeps_a_t_junction_stitched():
    """Plot Q's corners lie on plot P's left edge with no vertex there. When they move onto the
    settled boundary, P's edge must bend with them (46A opened a 100 m hairline this way)."""
    from shapely.geometry import LineString, Point
    A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])                       # settled
    P = Polygon([(10.3, 0), (20, 0), (20, 10), (10.3, 10)])                 # no vertices at y = 4, 7
    Q = Polygon([(10.3, 4), (14, 4), (14, 7), (10.3, 7)])                   # its west corners sit on P's edge
    P = P.difference(Q)                                                     # P wraps Q but keeps no vertices there? it does; so use the raw P below
    P_raw = Polygon([(10.3, 0), (20, 0), (20, 10), (10.3, 10)])
    out, _rep = topology.conform({"B": [({"poly_id": 1}, P_raw), ({"poly_id": 2}, Q)]}, {"A": A}, ["B"], anchors=())
    nP, nQ = out["B"][0][1], out["B"][1][1]
    assert nQ.distance(A) < 1e-9, "Q's west corners moved onto A"
    edge = LineString(nP.exterior.coords)
    for c in nQ.exterior.coords:
        if abs(c[0] - 10.0) < 1e-6:
            assert edge.distance(Point(c)) < 1e-6, "P's edge did not follow Q's moved corner %s" % (c,)


def test_node_parts_merges_near_duplicate_corners_and_inserts_t_junctions():
    P = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    Q = Polygon([(10.02, 3), (15, 3), (15, 6), (10.0, 6)])                   # (10.02,3) is 2 cm off P's corner line
    out, changed = topology.node_parts([P, Q])
    assert changed
    nP, nQ = out
    assert (10.0, 6.0) in {(round(x, 6), round(y, 6)) for x, y in nP.exterior.coords}, "T-junction vertex inserted into P"
    assert nP.is_valid and nQ.is_valid
