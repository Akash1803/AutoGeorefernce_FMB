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


def test_a_seam_between_two_parcels_goes_to_the_side_it_changes_least(monkeypatch):
    big = [("1", Polygon([(0, 0), (10, 0), (10, 20), (0, 20)])),
           ("2", Polygon([(10, 0), (20, 0), (20, 20), (10, 20)]))]
    small = [("1", Polygon([(0, 20.8), (20, 20.8), (20, 23.3), (0, 23.3)]))]
    done = {"BIG": big, "SMALL": small}
    monkeypatch.setattr(vc.neighbours, "printed", lambda v, s: ["SMALL"] if s == "BIG" else ["BIG"])
    vc.close_seams("test", done)
    b = unary_union([g for _l, g in done["BIG"]])
    s = unary_union([g for _l, g in done["SMALL"]])
    assert b.distance(s) < 1e-6, "the seam is closed"
    assert abs(s.area - 50.0) < 1e-6, "the small parcel keeps its drawing exactly"
    p1, p2 = dict(done["BIG"])["1"], dict(done["BIG"])["2"]
    assert abs(p1.bounds[2] - 10.0) < 1e-6 and abs(p2.bounds[0] - 10.0) < 1e-6, "the internal line runs straight on"
    assert abs(p1.bounds[3] - 20.8) < 1e-6 and abs(p2.bounds[3] - 20.8) < 1e-6


def test_a_seam_too_big_for_either_side_is_closed_down_its_middle(monkeypatch):
    a = [("1", Polygon([(0, 0), (20, 0), (20, 10), (0, 10)]))]
    b = [("1", Polygon([(0, 13), (20, 13), (20, 23), (0, 23)]))]          # a 3 m seam: 30 % of either
    done = {"A": a, "B": b}
    monkeypatch.setattr(vc.neighbours, "printed", lambda v, s: ["B"] if s == "A" else ["A"])
    vc.close_seams("test", done)
    ga, gb = done["A"][0][1], done["B"][0][1]
    assert ga.distance(gb) < 1e-6
    assert abs(ga.bounds[3] - 11.5) < 0.2 and abs(gb.bounds[1] - 11.5) < 0.2, "each side moves half"


def test_finish_coverage_makes_touching_lines_one_line_and_fills_pinholes():
    import shapely
    a = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    b = Polygon([(10.0004, 2), (20, 2), (20, 8), (10.0004, 8)])            # 0.4 mm off, no shared vertices
    c = Polygon([(0, 10.02), (10, 10.02), (10, 20), (0, 20)])              # 2 cm pin-gap to a
    done = {"A": [("", a)], "B": [("", b)], "C": [("", c)]}
    vc.finish_coverage(done)
    gs = [done[k][0][1] for k in ("A", "B", "C")]
    bad = shapely.coverage_invalid_edges(gs, gap_width=1.0)
    assert sum(x.length for x in bad if x is not None) < 1e-6, "a valid coverage: seams are exact"
    u = shapely.union_all(gs)
    assert u.geom_type == "Polygon" and len(u.interiors) == 0


def test_a_sliver_along_several_parcels_is_shared_out_not_given_to_one(monkeypatch):
    # three parcels in a row north of a long 0.4 m sliver, one long parcel south of it
    north = {"N%d" % i: [("", Polygon([(20 * i, 0.4), (20 * i + 20, 0.4), (20 * i + 20, 10), (20 * i, 10)]))]
             for i in range(3)}
    south = {"S": [("", Polygon([(0, -10), (60, -10), (60, 0), (0, 0)]))]}
    done = dict(north, **south)
    monkeypatch.setattr(vc.neighbours, "printed", lambda v, s: [])
    vc.close_seams("test", done)
    bodies = {s: unary_union([g for _l, g in p]) for s, p in done.items()}
    u = unary_union(list(bodies.values()))
    assert abs(u.area - (3 * 20 * 9.6 + 600 + 60 * 0.4)) < 0.05, "the sliver is gone"
    for s, b in bodies.items():
        assert b.bounds[2] - b.bounds[0] <= 60.0 + 1e-6
    for i in range(3):
        b = bodies["N%d" % i]
        assert b.bounds[0] >= 20 * i - 1e-6 and b.bounds[2] <= 20 * i + 20 + 1e-6, "no parcel grows a tail past its neighbours"


def test_tidy_outline_removes_a_needle_and_a_notch_but_keeps_corners_and_a_thin_strip():
    box = Polygon([(0, 0), (40, 0), (40, 20), (0, 20)])
    needle = Polygon([(40, 10), (46, 10.1), (40, 10.3)])
    notch = Polygon([(10, 19.4), (11, 19.4), (11, 20), (10, 20)])
    body = unary_union([box, needle]).difference(notch)
    out = vc.tidy_outline(body)
    assert abs(out.area - box.area) < 0.05
    assert out.bounds == (0.0, 0.0, 40.0, 20.0)
    strip = Polygon([(0, 0), (60, 0), (60, 2.6), (0, 2.6)])
    assert abs(vc.tidy_outline(strip).area - strip.area) < 1e-6, "a 2.6 m strip survives untouched"
    wedge = Polygon([(0, 0), (60, 0), (0, 8)])                       # an 8 degree tip
    assert abs(vc.tidy_outline(wedge).area - wedge.area) < 0.05, "a sharp but real tip survives"


def test_partition_tiles_the_outline_and_extends_internal_lines_straight():
    outline = Polygon([(0, 0), (20, 0), (20, 11), (0, 11)])
    plots = [("1", Polygon([(0, 0), (10.2, 0), (10.2, 10), (0, 10)])),      # overlaps plot 2 by 0.2 m
             ("2", Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]))]         # both stop 1 m short of the top
    out = vc.partition(plots, outline)
    u = unary_union([g for _l, g in out])
    assert abs(u.area - outline.area) < 1e-6
    assert out[0][1].intersection(out[1][1]).area < 1e-9
    p1, p2 = dict(out)["1"], dict(out)["2"]
    assert abs(p1.bounds[3] - 11) < 1e-9 and abs(p2.bounds[3] - 11) < 1e-9
    assert abs(p1.area - 110.0) < 2.5 and abs(p2.area - 110.0) < 2.5


def test_the_pdf_interior_is_refused_when_it_would_visibly_change_a_plot():
    cur = [("1", Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])), ("2", Polygon([(10, 0), (20, 0), (20, 10), (10, 10)]))]
    close = [("1", Polygon([(0, 0), (10.3, 0), (10.3, 10), (0, 10)])), ("2", Polygon([(10.3, 0), (20, 0), (20, 10), (10.3, 10)]))]
    far = [("1", Polygon([(0, 0), (13, 0), (13, 10), (0, 10)])), ("2", Polygon([(13, 0), (20, 0), (20, 10), (13, 10)]))]
    assert vc.interior_gate(close, cur)
    assert not vc.interior_gate(far, cur)
