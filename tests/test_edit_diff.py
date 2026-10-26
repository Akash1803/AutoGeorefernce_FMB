import geopandas as gpd
from shapely.geometry import Polygon

from autogeoref import edit_diff


def frame(geoms, surveys, plots):
    return gpd.GeoDataFrame({"survey_no": surveys, "plot_no": plots, "village_code": ["v"] * len(geoms)},
                            geometry=geoms, crs=32644)


A = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
B = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])


def test_an_edit_that_drags_a_neighbour_is_reported():
    before = frame([A, B], ["1", "2"], ["1", "1"])
    a2 = Polygon([(0, 0), (11, 0), (10, 10), (0, 10)])             # moved a corner shared with B
    b2 = Polygon([(11, 0), (20, 0), (20, 10), (10, 10)])           # ... and B moved with it
    rows = edit_diff.diff(before, frame([a2, b2], ["1", "2"], ["1", "1"]), expect={"1"})
    assert {r["survey"] for r in rows} == {"1", "2"}
    other = [r for r in rows if r["survey"] == "2"][0]
    assert not other["expected"] and other["of_which_shared_corners"] == 1


def test_untouched_plots_are_not_reported():
    before = frame([A, B], ["1", "2"], ["1", "1"])
    a2 = Polygon([(0, 0), (10, 0), (10, 10), (-1, 10)])            # its own corner only
    rows = edit_diff.diff(before, frame([a2, B], ["1", "2"], ["1", "1"]), expect={"1"})
    assert [r["survey"] for r in rows] == ["1"] and rows[0]["of_which_shared_corners"] == 0
