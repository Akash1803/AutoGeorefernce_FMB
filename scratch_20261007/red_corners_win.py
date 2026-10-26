"""Akash's OK (2026-10-07): where a red parcel and one of his meet at a corner and their PDFs disagree by less than
WIN_M, the red parcel's PDF corner wins. Corner to corner only: a shared corner lying more than OFF_M off the red
parcel's PDF line moves onto the red PDF corner within WIN_M of it (in every plot that has it); kink corners left
between a moved corner and the PDF line (within KINK_M of it) are dropped from every plot that has them.
Usage: python red_corners_win.py <out.gpkg>"""
import sys, os
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from autogeoref import ring_io, ring, fit, printed_numbers as pn, paths
from shape_check import fit_pdf
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
WIN_M, OFF_M, KINK_M = 1.5, 0.3, 3.0
V = "35_04_221"; K2 = 110540.0 * 108600.0
table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]; s = n.split("_parcels")[0]
    if "_parcels" not in n or s == "70A" or (s == "70B" and not n.endswith("_auto")):
        continue
    g = ring_io._read_source(L["source"])
    for pos in range(len(g)):
        rows.append({"layer": n, "survey": s, "red": n.endswith("_auto"), "pos": pos, "plot_no": str(g.plot_no.iat[pos]),
                     "geometry": g.geometry.iat[pos]})
A = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
U = A.to_crs(32644)
rings = [[tuple(p) for p in np.array(gm.exterior.coords)[:-1]] for gm in U.geometry]
orig_rings = [list(r) for r in rings]
owners = {}
for k, r in enumerate(rings):
    for p in r:
        owners.setdefault(p, set()).add(k)
moves, pdf_lines = {}, {}
for s in sorted(set(A.survey[A.red])):
    ks = [k for k in range(len(A)) if A.survey.iat[k] == s]
    cur = ring.largest(unary_union([U.geometry.iat[k] for k in ks]))
    sh = ring.load_sheet(V, s, table)
    iou, th, t = fit_pdf(sh.outline, cur)
    pdf = Polygon(fit.transform_points(sh.outline, th, t)); pdf_lines[s] = pdf.exterior
    corners = np.array(pdf.exterior.coords)[:-1]
    for k in ks:
        for v in rings[k]:
            shared_with_his = any(not A.red.iat[o] for o in owners[v])
            if not shared_with_his or v in moves:
                continue
            if pdf.exterior.distance(Point(v)) <= OFF_M:
                continue
            d = np.linalg.norm(corners - np.array(v), axis=1); c = int(np.argmin(d))
            if d[c] <= WIN_M:
                moves[v] = (tuple(corners[c]), s)
print("corner moves:", [(A.layer.iat[min(owners[v])].split("_parcels")[0], round(float(np.linalg.norm(np.subtract(t_, v))), 2), s) for v, (t_, s) in moves.items()])
# each move is tried alone and kept only if it opens no gap, makes no overlap and changes no plot by > MAX_M2
MAX_M2 = 5.0
def state_ok(rr):
    polys = [Polygon(r) for r in rr]
    if not all(p.is_valid for p in polys):
        return False
    t = shapely.STRtree(polys); a_, b_ = t.query(polys, predicate="intersects")
    for i, j in zip(a_, b_):
        if i < j and A.survey.iat[i] != A.survey.iat[j]:
            x = shapely.intersection(polys[i], polys[j])
            if sum(q.area for q in getattr(x, "geoms", [x]) if q.geom_type == "Polygon") > 1e-6:
                return False
    Ua = unary_union(polys)
    holes = [Polygon(r).difference(Ua).area for q in getattr(Ua, "geoms", [Ua]) for r in q.interiors]
    return sum(1 for h in holes if h > 1e-6) <= BASE_HOLES
_U0 = unary_union([Polygon(r) for r in rings])
BASE_HOLES = sum(1 for q in getattr(_U0, "geoms", [_U0]) for r in q.interiors if Polygon(r).difference(_U0).area > 1e-6)
kept = {}
for v, (t_, s) in moves.items():
    trial = [list(r) for r in rings]
    for k in owners[v]:
        trial[k] = [t_ if p == v else p for p in trial[k]]
    # the kink corner beside it (shared with his, off the red PDF line, within KINK_M) goes too
    for w, ks in owners.items():
        if w != v and owners[w] & owners[v] and np.linalg.norm(np.subtract(w, t_)) <= KINK_M                 and pdf_lines[s].distance(Point(w)) > OFF_M and any(A.red.iat[k] for k in ks) and not all(A.red.iat[k] for k in ks):
            for k in ks:
                trial[k] = [p for p in trial[k] if p != w]
    da = max(abs(Polygon(trial[k]).area - Polygon(rings[k]).area) for k in range(len(rings)) if trial[k] != rings[k]) if any(trial[k] != rings[k] for k in range(len(rings))) else 0
    ok = da <= MAX_M2 and state_ok(trial)
    print("  move %.2f m at %s corner: largest plot change %.2f m2 -> %s" % (float(np.linalg.norm(np.subtract(t_, v))), s, da, "kept" if ok else "skipped"))
    if ok:
        rings = trial; kept[v] = (t_, s)
moves = kept
# kink corners near a moved corner, off the red PDF line, shared by the red parcel and his: dropped everywhere
moved_pts = []
drops = set()
for v, ks in owners.items():
    if v in moves or not moved_pts:
        continue
    reds = {A.survey.iat[k] for k in ks if A.red.iat[k]}
    if not reds or all(A.red.iat[k] for k in ks):
        continue
    near = min(np.linalg.norm(m - np.array(v)) for m in moved_pts)
    if near <= KINK_M and all(pdf_lines[s].distance(Point(v)) > OFF_M for s in reds):
        drops.add(v)
for k in range(len(rings)):
    r = [p for p in rings[k] if p not in drops]
    rings[k] = [p for i, p in enumerate(r) if p != r[i - 1]]
print("kink corners dropped:", len(drops))
new = []
for k in range(len(rings)):
    p = Polygon(rings[k])
    new.append(p if p.is_valid else U.geometry.iat[k])
    if not p.is_valid: print("  kept original (invalid after edit):", A.layer.iat[k], A.plot_no.iat[k])
# back to EPSG:4326; untouched corners keep their exact original coordinates
ll = list(gpd.GeoSeries(new, crs=32644).to_crs(4326))
from scipy.spatial import cKDTree
opts = np.vstack([shapely.get_coordinates(g_) for g_ in A.geometry]); kd = cKDTree(opts)
def back(gm):
    def f(c):
        d, k = kd.query(c); c2 = c.copy(); near = d <= 1e-8; c2[near] = opts[k[near]]; return c2
    return shapely.transform(gm, f)
B = A.copy(); B["geometry"] = [A.geometry.iat[k] if rings[k] == orig_rings[k] else back(ll[k]) for k in range(len(A))]
ch = [(B.layer.iat[k], B.plot_no.iat[k], round((gpd.GeoSeries([B.geometry.iat[k]], crs=4326).to_crs(32644).area.iat[0] - U.geometry.iat[k].area), 2))
      for k in range(len(B)) if rings[k] != orig_rings[k]]
print("plots changed (area change m2):", ch)
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); ov = []
for i, j in zip(a_, b_):
    if i < j and B.survey.iat[i] != B.survey.iat[j]:
        x = shapely.intersection(g[i], g[j]); ar = sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon")
        if ar > 0: ov.append((B.survey.iat[i], B.survey.iat[j], ar * K2))
Ua = unary_union(g); holes = [round(Polygon(r).difference(Ua).area * K2, 3) for p in getattr(Ua, "geoms", [Ua]) for r in p.interiors]
print("overlaps", ov, "| holes", holes, "| invalid", int((~B.is_valid).sum()))
B.drop(columns=["red"]).to_file(OUT, layer="plan", driver="GPKG")
