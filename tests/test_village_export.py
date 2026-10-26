import json
import zipfile

import geopandas as gpd
from shapely.geometry import Polygon

from autogeoref import paths, village_export, visible


def test_export_reads_only_what_qgis_shows(tmp_path, monkeypatch):
    v = "35_04_074"
    vd = tmp_path / "FMB_Vector" / v
    vd.mkdir(parents=True)
    his = Polygon([(391000, 1411000), (391030, 1411000), (391030, 1411020), (391000, 1411020)])
    hidden = Polygon([(391500, 1411000), (391530, 1411000), (391530, 1411020), (391500, 1411020)])
    f = vd / "565_parcels_modified.gpkg"
    gpd.GeoDataFrame({"plot_no": ["2A"]}, geometry=[his], crs=32644).to_file(f, layer="565_parcels_modified", driver="GPKG")
    gpd.GeoDataFrame({"plot_no": ["2A"]}, geometry=[hidden], crs=32644).to_file(f, layer="parcels", driver="GPKG")
    sheet = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {"plot_no": "2A", "poly_id": 0},
             "geometry": {"type": "Polygon", "coordinates": [[[0, 0], [30, 0], [30, 20], [0, 20], [0, 0]]]}}]}
    (vd / "565_parcels.geojson").write_text(json.dumps(sheet) + " " * 300)
    proj = tmp_path / "p.qgz"
    with zipfile.ZipFile(proj, "w") as z:
        z.writestr("p.qgs", "<qgis><maplayer><datasource>FMB_Vector/35_04_074/565_parcels_modified.gpkg"
                            "|layername=565_parcels_modified</datasource><layername>565</layername></maplayer></qgis>")
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    monkeypatch.setattr(paths, "QGIS_PROJECT", proj)
    visible._project_sources.cache_clear()
    got = village_export.build(v, "Thirukatchur", tmp_path / "out", workbook=tmp_path / "none.xlsx")
    plots = gpd.read_file(got["plots"]).to_crs(32644)
    assert len(plots) == 1 and plots.geometry.iloc[0].centroid.distance(his.centroid) < 0.01
    assert list(plots.columns[:4]) == ["kide", "survey_no", "subdiv_no", "area_acre"]
    assert plots.kide.iloc[0] == "565/2A" and plots.iou_vs_pdf.iloc[0] > 0.99
