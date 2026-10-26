"""Read-only: per survey, PDF title-block area (read twice), converted vector area, placed geometry area
(measured twice: shapely in UTM 44N, QGIS ellipsoidal), plus plot-level vector vs geometry."""
import sys, os, json, glob
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, geopandas as gpd
from shapely.ops import unary_union
from autogeoref import printed_numbers as pn, qgis_bridge
PR = r"D:\Projects\Tambaram_Chengalpattu"; V = "35_04_221"
# second PDF reading: by eye from the cropped title blocks (ares; hectares x 100)
EYE = {"264": 79.00, "71B": 50.00, "71A": 3.00, "72": 39.00, "74": 53.30, "66A": 121.50, "66B": 20.50, "67A": 67.50,
       "67B": 33.00, "68A": 18.00, "68B": 18.00, "69": 7.50, "70A": 25.50, "70B": 102.50, "265": 26.00, "266": 39.50, "267": 52.00}
AUTO = ["264", "71B", "71A", "72", "74"]
rows, plots = [], []
for s in EYE:
    auto = s in AUTO
    f = os.path.join(PR, "FMB_Vector", V, s + "_parcels_auto.gpkg") if auto else os.path.join(PR, s + "_parcels_modified.geojson")
    G = gpd.read_file(f, layer="parcels") if auto else gpd.read_file(f)
    vec = json.load(open(os.path.join(PR, "FMB_Vector", V, s + "_parcels.geojson"), encoding="utf-8"))
    vec_area = sum(ft["properties"]["area_sqm"] for ft in vec["features"])
    geo = G.to_crs(32644)
    rows.append({"survey": s, "who": "auto (red)" if auto else "yours", "pdf_read_1": round(pn.read_header_area(os.path.join(PR, "FMB_Sketches", V, s + ".pdf")) or 0, 1),
                 "pdf_read_2": round(EYE[s] * 100, 1), "vector_m2": round(vec_area, 1),
                 "geometry_m2": round(unary_union(geo.geometry.values).area, 1), "plots": len(G), "vector_plots": len(vec["features"])})
    vp = {ft["properties"]["poly_id"]: ft["properties"]["area_sqm"] for ft in vec["features"]}
    for _i, r in geo.iterrows():
        plots.append({"survey": s, "who": rows[-1]["who"], "poly_id": r.get("poly_id"), "plot_no": r.get("plot_no"),
                      "vector_m2": round(vp.get(r.get("poly_id"), float("nan")), 1), "geometry_m2": round(r.geometry.area, 1)})
D = pd.DataFrame(rows)
code = """
import json
from qgis.core import QgsProject, QgsDistanceArea, QgsCoordinateReferenceSystem
P = QgsProject.instance(); d = QgsDistanceArea(); d.setEllipsoid("WGS84")
sub = P.layerTreeRoot().findGroup("Manual_Georeferenced").findGroup("Paranur_Village")
out = {}
for n in sub.findLayers():
    L = n.layer()
    if L is None or "_parcels" not in L.name(): continue
    d.setSourceCrs(L.crs(), P.transformContext())
    s = L.name().split("_parcels")[0]
    out[s] = sum(d.measureArea(f.geometry()) for f in L.getFeatures())
print(json.dumps(out))
"""
q = qgis_bridge.json_result(code)
D["geometry_qgis_m2"] = D.survey.map(lambda s: round(q.get(s, float("nan")), 1))
D["pdf_reads_agree"] = (D.pdf_read_1 - D.pdf_read_2).abs() < 0.6
D["vector_vs_pdf_%"] = (100 * (D.vector_m2 / D.pdf_read_2 - 1)).round(1)
D["geometry_vs_vector_%"] = (100 * (D.geometry_m2 / D.vector_m2 - 1)).round(1)
D["geometry_vs_pdf_%"] = (100 * (D.geometry_m2 / D.pdf_read_2 - 1)).round(1)
D["two_geometry_measures_%"] = (100 * (D.geometry_qgis_m2 / D.geometry_m2 - 1)).round(2)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 30)
print(D.to_string(index=False))
P = pd.DataFrame(plots)
P["diff_%"] = (100 * (P.geometry_m2 / P.vector_m2 - 1)).round(1)
print("\nplots, vector vs geometry (auto parcels):")
print(P[P.who == "auto (red)"].to_string(index=False))
print("\nplots beyond +-5 %% (yours):")
print(P[(P.who == "yours") & (P["diff_%"].abs() > 5)].to_string(index=False))
out = os.path.join(sorted(glob.glob(os.path.join(PR, "Auto_Georef", V, "batch_*")))[-1], "area_check.xlsx")
with pd.ExcelWriter(out) as xw:
    D.to_excel(xw, sheet_name="surveys", index=False); P.to_excel(xw, sheet_name="plots", index=False)
print(out)
