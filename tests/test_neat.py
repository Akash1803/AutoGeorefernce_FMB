import geopandas as gpd
import shapely
from shapely.geometry import Polygon

from autogeoref import neat


def box(x0, y0, x1, y1):
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def frame(geoms, surveys):
    return gpd.GeoDataFrame({"survey_no": surveys}, geometry=geoms, crs=32644)


def on_old_lines(out, before, tol=1e-3):
    """Every vertex of the result lies on a line that was already drawn: no new boundaries."""
    old = shapely.union_all([g.boundary for g in before])
    for g in out.geometry:
        for x, y in shapely.get_coordinates(g):
            if old.distance(shapely.Point(x, y)) > tol:
                return False
    return True


def test_an_overlap_is_closed_on_an_existing_line():
    before = [box(0, 0, 10, 10), box(9, 0, 20, 10)]
    out, _ = neat.clean(frame(before, ["1", "2"]), {"1": {"2"}, "2": {"1"}})
    a, b = out.geometry
    assert a.intersection(b).area < 1e-3
    assert abs(a.area + b.area - 200) < 0.01                          # no ground lost or added
    assert shapely.coverage_is_valid(list(out.geometry))
    assert on_old_lines(out, before)


def test_a_splay_between_pdf_neighbours_is_closed_by_moving_a_corner():
    a = box(0, 0, 10, 10)
    b = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])               # touches a at one corner only
    out, _ = neat.clean(frame([a, b], ["1", "2"]), {"1": {"2"}, "2": {"1"}})
    u = shapely.union_all(list(out.geometry))
    assert abs(u.area - (a.area + b.area + 15)) < 0.5                  # the 15 m2 wedge is closed
    assert out.geometry.iloc[0].intersection(out.geometry.iloc[1]).area < 1e-3
    assert all(len(shapely.get_coordinates(g)) <= 6 for g in out.geometry)   # still simple shapes


def test_a_real_lane_between_non_neighbours_is_kept():
    a, b = box(0, 0, 10, 10), box(14, 0, 24, 10)                      # a 4 m lane, not PDF neighbours
    out, _ = neat.clean(frame([a, b], ["1", "2"]), {})
    assert out.geometry.iloc[0].equals(a) and out.geometry.iloc[1].equals(b)


def test_a_hairline_between_any_two_plots_is_closed():
    before = [box(0, 0, 10, 10), box(10.4, 0, 20, 10)]                # 0.4 m double line
    out, _ = neat.clean(frame(before, ["1", "2"]), {})
    assert out.geometry.iloc[0].distance(out.geometry.iloc[1]) < 1e-3
    assert shapely.coverage_is_valid(list(out.geometry))
    assert on_old_lines(out, before)


def test_plots_far_from_any_seam_keep_their_exact_vertices():
    a, b, c = box(0, 0, 10, 10), box(9, 0, 20, 10), box(100, 100, 110, 110)
    out, _ = neat.clean(frame([a, b, c], ["1", "2", "3"]), {"1": {"2"}, "2": {"1"}})
    assert out.geometry.iloc[2].equals(c)


def test_a_thin_strip_is_never_collapsed_by_its_neighbours():
    left, strip, right = box(0, 0, 10, 20), box(10.3, 0, 12.8, 20), box(13.1, 0, 23, 20)   # 2.5 m strip
    out, _ = neat.clean(frame([left, strip, right], ["1", "2", "3"]),
                        {"1": {"2"}, "2": {"1", "3"}, "3": {"2"}})
    assert out.geometry.iloc[1].area > 0.9 * strip.area
    assert shapely.coverage_is_valid(list(out.geometry))


def test_no_new_small_pieces_are_created():
    a1, a2 = box(0, 0, 10, 5), box(0, 5, 10, 10)
    b = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])
    out, _ = neat.clean(frame([a1, a2, b], ["1", "1", "2"]), {"1": {"2"}, "2": {"1"}})
    assert all(g.geom_type == "Polygon" for g in out.geometry)
    assert shapely.coverage_is_valid(list(out.geometry))


def test_a_vertex_never_snaps_across_another_parcel():
    big = box(0, 0, 30, 30)
    strip = box(30, 0, 31.5, 30)                                   # a 1.5 m strip between them
    other = Polygon([(31.5, 0), (40, 0), (40, 30), (33, 30)])       # splays away from the strip
    out, _ = neat.clean(frame([big, strip, other], ["1", "2", "3"]),
                        {"1": {"2", "3"}, "2": {"1", "3"}, "3": {"1", "2"}})
    assert out.geometry.iloc[1].area > 0.95 * strip.area             # the strip keeps its ground
    assert out.geometry.iloc[0].intersection(out.geometry.iloc[1]).area < 1e-3


def test_an_internal_line_is_extended_along_itself_not_swung():
    # survey 1 = two plots split by a slanted line; survey 2 sits 2 m away on a straight edge
    p1 = Polygon([(0, 0), (6, 0), (4, 10), (0, 10)])
    p2 = Polygon([(6, 0), (10, 0), (10, 10), (4, 10)])
    other = box(0, 12, 10, 20)                                          # 2 m gap above survey 1
    out, _ = neat.clean(frame([other, p1, p2], ["2", "1", "1"]), {"1": {"2"}, "2": {"1"}})
    a, b = out.geometry.iloc[1], out.geometry.iloc[2]
    seam = a.intersection(b)
    xs = sorted(x for x, _y in shapely.get_coordinates(seam))
    # the slanted line keeps its slope: at y = 12 it reaches x = 6 - 12 / 5 = 3.6
    top = max(shapely.get_coordinates(seam), key=lambda c: c[1])
    assert abs(top[1] - 12) < 0.01 and abs(top[0] - 3.6) < 0.05


def test_a_corner_two_surveys_already_share_is_never_dragged():
    # surveys 1 and 3 meet at (10, 10); survey 2 sits 3 m away: closing it must not move that corner
    s1 = box(0, 0, 10, 10)
    s3 = box(0, 10, 10, 20)
    s2 = Polygon([(13, 0), (23, 0), (23, 20), (13, 20)])
    before = [s2, s1, s3]
    out, _ = neat.clean(frame(before, ["2", "1", "3"]), {"1": {"2"}, "2": {"1", "3"}, "3": {"2"}})
    corner = shapely.Point(10, 10)
    assert out.geometry.iloc[1].boundary.distance(corner) < 1e-6      # still on survey 1
    assert out.geometry.iloc[2].boundary.distance(corner) < 1e-6      # still on survey 3
    assert abs(out.geometry.iloc[1].area - s1.area) < 1e-6 or out.geometry.iloc[1].area > s1.area
