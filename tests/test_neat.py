import geopandas as gpd
import shapely
from shapely.geometry import Polygon

from autogeoref import neat


def box(x0, y0, x1, y1):
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def frame(geoms, surveys):
    return gpd.GeoDataFrame({"survey_no": surveys}, geometry=geoms, crs=32644)


def test_an_overlap_is_split_along_its_middle_into_one_line():
    g = frame([box(0, 0, 10, 10), box(9, 0, 20, 10)], ["1", "2"])
    out, rep = neat.clean(g, {"1": {"2"}, "2": {"1"}})
    a, b = out.geometry
    assert a.intersection(b).area < 1e-3
    assert abs(a.area - 95) < 0.2 and abs(b.area - 105) < 0.2        # each yields half the strip
    assert shapely.coverage_is_valid(list(out.geometry))


def test_a_splay_between_pdf_neighbours_is_closed():
    a = box(0, 0, 10, 10)
    b = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])               # touches a at one corner only
    out, _ = neat.clean(frame([a, b], ["1", "2"]), {"1": {"2"}, "2": {"1"}})
    u = shapely.union_all(list(out.geometry))
    assert abs(u.area - (a.area + b.area + 15)) < 0.5                  # the 15 m2 wedge is filled
    assert out.geometry.iloc[0].intersection(out.geometry.iloc[1]).area < 1e-3


def test_a_real_lane_between_non_neighbours_is_kept():
    a, b = box(0, 0, 10, 10), box(14, 0, 24, 10)                      # a 4 m lane, not PDF neighbours
    out, _ = neat.clean(frame([a, b], ["1", "2"]), {})
    assert abs(out.geometry.iloc[0].area - 100) < 1e-6 and abs(out.geometry.iloc[1].area - 100) < 1e-6


def test_a_hairline_between_any_two_plots_is_closed():
    a, b = box(0, 0, 10, 10), box(10.4, 0, 20, 10)                    # 0.4 m double line
    out, _ = neat.clean(frame([a, b], ["1", "2"]), {})
    assert out.geometry.iloc[0].distance(out.geometry.iloc[1]) < 1e-3
    assert shapely.coverage_is_valid(list(out.geometry))


def test_plots_far_from_any_seam_keep_their_exact_vertices():
    a, b, c = box(0, 0, 10, 10), box(9, 0, 20, 10), box(100, 100, 110, 110)
    out, _ = neat.clean(frame([a, b, c], ["1", "2", "3"]), {"1": {"2"}, "2": {"1"}})
    assert out.geometry.iloc[2].equals(c)


def test_a_wedge_where_three_plots_meet_is_shared_among_them():
    a1, a2 = box(0, 0, 10, 5), box(0, 5, 10, 10)                      # survey 1 in two plots
    b = Polygon([(10, 10), (13, 0), (23, 0), (20, 10)])               # survey 2 splays away from both
    out, _ = neat.clean(frame([a1, a2, b], ["1", "1", "2"]), {"1": {"2"}, "2": {"1"}})
    u = shapely.union_all(list(out.geometry))
    assert abs(u.area - (a1.area + a2.area + b.area + 15)) < 0.5
    assert out.geometry.iloc[0].area > 50.5 and out.geometry.iloc[1].area > 50.5
    assert shapely.coverage_is_valid(list(out.geometry))
