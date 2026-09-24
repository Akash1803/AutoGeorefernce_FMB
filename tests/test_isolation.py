"""The automation and Akash's hand files never share a file (2026-09-24 recovery).

On 23 Sep the engine wrote its placements as <s>_parcels_modified.gpkg in his folder, scripts
added hidden layers to his files, and approved parcels existed twice. These tests pin the rules:
the tool writes only into its own _auto folder, never into a file his QGIS project shows, and a
parcel his project shows is always an anchor, whoever first drew it.
"""
import zipfile

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Polygon

from autogeoref import anchors, paths, pose_solver, visible

SQ = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])


def _gpkg(path, layer="parcels"):
    gpd.GeoDataFrame({"plot_no": ["1"], "poly_id": [0]}, geometry=[SQ], crs=32644).to_file(
        path, layer=layer, driver="GPKG")


def _project(tmp_path, monkeypatch, village, shown):
    vd = tmp_path / "FMB_Vector" / village
    vd.mkdir(parents=True, exist_ok=True)
    proj = tmp_path / "CUMTA" / "p.qgz"
    proj.parent.mkdir(exist_ok=True)
    xml = "<qgis>" + "".join(
        "<maplayer><datasource>../FMB_Vector/%s/%s|layername=%s</datasource><layername>%s</layername></maplayer>"
        % (village, f, lay, f) for f, lay in shown) + "</qgis>"
    with zipfile.ZipFile(proj, "w") as z:
        z.writestr("p.qgs", xml)
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    monkeypatch.setattr(paths, "QGIS_PROJECT", proj)
    visible._project_sources.cache_clear()
    return vd


def test_the_tool_writes_only_into_its_own_folder(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    out = paths.output_path("35_04_074", "53")
    assert out.parent == paths.vector_dir("35_04_074") / "_auto"
    assert out.name == "53_parcels_auto.gpkg"
    assert out not in paths.manual_files("35_04_074", "53")


def test_a_file_his_project_shows_is_never_written_by_the_tool(tmp_path, monkeypatch):
    vd = _project(tmp_path, monkeypatch, "35_04_074", [("53_parcels_auto.gpkg", "parcels")])
    _gpkg(vd / "53_parcels_auto.gpkg")
    with pytest.raises(RuntimeError):
        visible.refuse_tool_write(vd / "53_parcels_auto.gpkg")     # approved and shown: his now
    visible.refuse_tool_write(paths.output_path("35_04_074", "53"))  # the tool's own copy: fine


def test_a_parcel_his_project_shows_is_always_an_anchor(tmp_path, monkeypatch):
    vd = _project(tmp_path, monkeypatch, "35_04_074", [("53_parcels_auto.gpkg", "parcels")])
    f = vd / "53_parcels_auto.gpkg"
    _gpkg(f)
    monkeypatch.setattr(paths, "surveys_with_sheets", lambda v: ["53"])
    monkeypatch.setattr(anchors, "pose_from_points", lambda v, s: (0.0, np.zeros(2), 0.1, 1.0, 1.0, 3))
    monkeypatch.setattr(anchors, "_area_scale", lambda v, s, f: 1.0)
    fp = anchors.fingerprint(f)
    got = anchors.load_anchors("35_04_074", {fp})                    # the tool drew it first
    assert "53" in got and got["53"].file == f


def test_pose_solver_never_applies_into_his_files():
    with pytest.raises(RuntimeError):
        pose_solver.apply("35_04_074")
