import zipfile
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import LineString, Polygon

from autogeoref import paths, visible


def _gpkg(path, layers):
    for name, geom in layers:
        gpd.GeoDataFrame({"plot_no": ["1"]}, geometry=[geom], crs=32644).to_file(path, layer=name, driver="GPKG")


SQ = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
MOVED = Polygon([(5, 0), (15, 0), (15, 10), (5, 10)])


def test_his_named_layer_wins_over_a_hidden_tool_copy(tmp_path):
    f = tmp_path / "565_parcels_modified.gpkg"
    _gpkg(f, [("565_parcels_modified", SQ), ("parcels", MOVED)])
    assert visible.hand_layer(f) == "565_parcels_modified"
    assert visible.has_hidden_copy(f)
    assert visible.read_hand(f).geometry.iloc[0].equals(SQ)


def test_edges_listed_first_are_never_read_as_the_parcel(tmp_path):
    f = tmp_path / "54_parcels_modified.gpkg"
    _gpkg(f, [("edges", LineString([(0, 0), (10, 0)])), ("parcels", SQ)])
    assert visible.hand_layer(f) == "parcels"
    assert visible.read_hand(f).geometry.iloc[0].geom_type == "Polygon"


def test_the_tool_refuses_to_write_parcels_into_a_hand_file(tmp_path):
    hand = tmp_path / "609_parcels_modified.gpkg"
    _gpkg(hand, [("609_parcels_modified", SQ)])
    with pytest.raises(RuntimeError):
        visible.refuse_tool_write(hand)
    tool = tmp_path / "39A_parcels_modified.gpkg"
    _gpkg(tool, [("parcels", SQ), ("edges", LineString([(0, 0), (10, 0)]))])
    visible.refuse_tool_write(tool)                      # its own file: fine
    visible.refuse_tool_write(tmp_path / "new_parcels_modified.gpkg")


def test_shown_reads_the_file_and_layer_from_the_project(tmp_path, monkeypatch):
    vd = tmp_path / "FMB_Vector" / "35_04_074"
    vd.mkdir(parents=True)
    _gpkg(vd / "53_parcels_auto.gpkg", [("parcels", SQ)])
    _gpkg(vd / "53_parcels_modified.gpkg", [("parcels", MOVED)])
    _gpkg(vd / "565_parcels_modified.gpkg", [("565_parcels_modified", SQ), ("parcels", MOVED)])
    proj = tmp_path / "CUMTA" / "p.qgz"
    proj.parent.mkdir()
    xml = "<qgis>" + "".join(
        "<maplayer><datasource>../FMB_Vector/35_04_074/%s|layername=%s</datasource><layername>%s</layername></maplayer>"
        % (f, lay, name) for f, lay, name in [("53_parcels_auto.gpkg", "parcels", "53_parcels_approved"),
                                              ("565_parcels_modified.gpkg", "565_parcels_modified", "565")]) + "</qgis>"
    with zipfile.ZipFile(proj, "w") as z:
        z.writestr("p.qgs", xml)
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    monkeypatch.setattr(paths, "QGIS_PROJECT", proj)
    sh = visible.shown("35_04_074")
    assert sh["53"][0].name == "53_parcels_auto.gpkg" and sh["53"][1] == "parcels"
    assert sh["565"][1] == "565_parcels_modified"
    assert paths.manual_files("35_04_074", "53")[0].name == "53_parcels_auto.gpkg"
