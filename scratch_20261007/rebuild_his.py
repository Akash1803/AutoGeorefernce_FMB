"""Step A of the clean rebuild (plan only; nothing of his is written here):
* 267, 266, 68B, 69 (accepted PDF versions): back to the exact PDF drawing (rigid pose fitted to the current shape);
  where two of these drawings disagree the later one gives way along the earlier one's line;
* 66A, 66B, 67A, 67B, 68A, 265: his originals from before today's edits, giving way to those four only by an exact
  cut along their lines (and taking gaps <= 3 m wide) - no vertex is moved;
* every seam then carries the same vertices on both sides (collinear points inserted).
Usage: python rebuild_his.py <out.gpkg>"""
import sys, os
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from autogeoref import ring_io, ring, fit, ring_topo, printed_numbers as pn, paths
from shape_check import fit_pdf
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
V = "35_04_221"; PR = paths.PROJECT
BACKUP = PR / "_logs" / "paranur_accept_pdf4_20261007_1640"
ACCEPT = ["267", "266", "68B", "69"]                       # fixed first, in this order
GIVE = ["67A", "67B", "68A", "66A", "66B", "265"]
table = pn.village_table(V, str(PR / "Auto_Georef" / V / "printed_numbers.json"))
src = {L["name"]: L["source"] for L in ring_io.village_layers("Paranur")}
rows, fixed = [], []
for s in ACCEPT:
    cur = ring_io._read_source(src[s + "_parcels_modified"])
    assert len(cur) == 1
    sh = ring.load_sheet(V, s, table)
    target = ring.largest(unary_union(cur.to_crs(32644).geometry.values))
    iou, th, t = fit_pdf(sh.outline, target)
    exact = gpd.GeoSeries([fit.apply_pose(sh.plots[0][1] if len(sh.plots) == 1 else sh.polygon, th, t)], crs=32644).to_crs(4326).iat[0]
    g = ring_topo.exact_seams([exact], fixed, 3.0)[0] if fixed else exact
    print("%s: PDF fit IoU %.4f vs current; gave way %.2f m2 to earlier PDF parcels" % (s, iou, (exact.area - g.area) * 110540 * 108600))
    fixed.append(g)
    rows.append({"layer": s + "_parcels_modified", "survey": s, "pos": 0, "plot_no": str(cur.plot_no.iat[0]), "geometry": g})
for s in GIVE:
    orig = gpd.read_file(BACKUP / ("%s_parcels_modified.geojson" % s))
    plots = list(orig.geometry)
    new = ring_topo.exact_seams(plots, fixed, 3.0)
    for pos, g in enumerate(new):
        rows.append({"layer": s + "_parcels_modified", "survey": s, "pos": pos, "plot_no": str(orig.plot_no.iat[pos]), "geometry": g})
B = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
# same vertices on both sides of every seam
geoms = list(B.geometry)
for i in range(len(geoms)):
    others = [geoms[j] for j in range(len(geoms)) if j != i and geoms[j].distance(geoms[i]) < 1e-8]
    geoms[i] = ring_topo.insert_nodes(geoms[i], others)
B["geometry"] = geoms
K2 = 110540.0 * 108600.0
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects")
ov = []
for i, j in zip(a_, b_):
    if i < j and B.survey.iat[i] != B.survey.iat[j]:
        x = shapely.intersection(g[i], g[j]); ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
        if ar > 0: ov.append((B.survey.iat[i], B.survey.iat[j], ar * K2))
U = unary_union(g); holes = [round(Polygon(r).difference(U).area * K2, 3) for p in getattr(U, "geoms", [U]) for r in p.interiors]
print("his set: overlaps", ov, "| holes", holes, "| invalid", int((~B.is_valid).sum()))
B.to_file(OUT, layer="plan", driver="GPKG")
print("written", OUT)
