from shapely.geometry import Polygon

from autogeoref import fabric_check as fc


def sq(x0, y0, w=10.0, h=10.0):
    return Polygon([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)])


def test_overlaps_are_measured_between_surveys():
    got = fc.overlaps({"1": sq(0, 0), "2": sq(9, 0), "3": sq(30, 0)})
    assert [(a, b) for a, b, _ in got] == [("1", "2")]
    assert abs(got[0][2] - 10.0) < 1e-6


def test_a_shared_edge_is_closed_and_a_splay_is_an_open_gap():
    a = sq(0, 0)
    closed = sq(10, 0)
    # touches a at one corner only and opens a wedge along the shared side
    splay = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])
    assert fc.seam_gap(a, closed) < 0.3
    assert fc.seam_gap(a, splay) > 5.0            # min distance says 0.00 m; the wedge is real


def test_gaps_only_between_pdf_neighbours():
    bodies = {"1": sq(0, 0), "2": Polygon([(10, 10), (13, 0), (23, 0), (20, 10)]), "3": sq(40, 0)}
    nbrs = {"1": {"2"}, "2": {"1"}}
    got = fc.gaps(bodies, nbrs)
    assert [(a, b) for a, b, _m2, _note in got] == [("1", "2")]


def test_rigid_fit_scores_the_shape_against_the_pdf():
    drawing = Polygon([(0, 0), (30, 0), (30, 10), (10, 20), (0, 20)])
    import shapely.affinity as af
    placed = af.translate(af.rotate(drawing, 37, origin=(0, 0)), 391000, 1411000)
    iou, pose = fc.shape_iou(drawing, placed)
    assert iou > 0.99
    stretched = af.scale(placed, 1.15, 1.0)
    assert fc.shape_iou(drawing, stretched)[0] < 0.95


def test_ground_covered_by_a_third_parcel_is_not_a_gap():
    a = sq(0, 0)
    b = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])
    wedge_filler = Polygon([(10, 0), (13, 0), (10, 10)])            # a third survey sits in the wedge
    got = fc.gaps({"1": a, "2": b, "3": wedge_filler}, {"1": {"2"}, "2": {"1"}})
    assert got == []
