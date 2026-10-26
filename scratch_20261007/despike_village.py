"""Remove spikes, hooks and tiny steps from every Paranur parcel, seam by seam so neighbours stay identical.

All parcel boundaries (his + red, except his old 70A/70B which the red PDF versions are about to replace)
are noded into arcs between junctions; each arc is cleaned once (needle spikes removed, steps shorter than
STEP_M collapsed, then Douglas-Peucker at TOL_M, junctions kept), and the parcels are rebuilt from the
cleaned arcs. Untouched vertices keep their exact original coordinates. Writes a plan to OUT only.
Usage: python despike_village.py <out.gpkg>"""
import sys, os, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union, linemerge, polygonize
from scipy.spatial import cKDTree
from autogeoref import ring_io
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
SKIP_HIS = {"70A", "70B"}           # his old versions, pending replacement by the red PDF versions
SPIKE_DEG, SPIKE_EDGE_M = 25.0, 6.0
STEP_M = 1.0
TOL_M = 0.3

rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]
    if "_parcels" not in n:
        continue
    s = n.split("_parcels")[0]
    red = n.endswith("_parcels_auto")
    if s in SKIP_HIS and not red:
        continue
    g = ring_io._read_source(L["source"])
    for pos in range(len(g)):
        if g.geometry.iat[pos] is None or g.geometry.iat[pos].is_empty:
            continue
        rows.append({"layer": n, "survey": s, "pos": pos, "plot_no": str(g.plot_no.iat[pos]) if "plot_no" in g else "",
                     "geometry": g.geometry.iat[pos]})
A = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
U = A.to_crs(32644)
orig = list(U.geometry)

# 1. noded arcs between junctions
lines = unary_union([p.exterior for g_ in orig for p in getattr(g_, "geoms", [g_])])
arcs = list(getattr(linemerge(lines), "geoms", [linemerge(lines)]))

def clean_arc(c):
    c = [tuple(p) for p in c]
    closed = math.dist(c[0], c[-1]) < 1e-9
    changed = True
    while changed and len(c) > 2:
        changed = False
        for i in range(1, len(c) - 1):                    # endpoints (junctions) never move
            a, b, d = np.array(c[i - 1]), np.array(c[i]), np.array(c[i + 1])
            e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
            if e1 < 1e-6 or e2 < 1e-6:
                del c[i]; changed = True; break
            ang = math.degrees(math.acos(max(-1, min(1, np.dot(a - b, d - b) / (e1 * e2)))))
            if ang < SPIKE_DEG and min(e1, e2) < SPIKE_EDGE_M:
                del c[i]; changed = True; break         # a needle: the tip goes
        if changed:
            continue
        for i in range(1, len(c) - 2):                    # a step: short edge between two turns
            b, d = np.array(c[i]), np.array(c[i + 1])
            if np.linalg.norm(d - b) < STEP_M:
                a, e = np.array(c[i - 1]), np.array(c[i + 2])
                if np.linalg.norm(b - a) > STEP_M and np.linalg.norm(e - d) > STEP_M:
                    mid = tuple((b + d) / 2)
                    c[i:i + 2] = [mid]; changed = True; break
    ls = LineString(c)
    s = ls.simplify(TOL_M, preserve_topology=False)
    out = list(s.coords)
    if closed and len(out) < 4:
        out = c
    return out

new_arcs = []
for a in arcs:
    c = list(a.coords)
    new_arcs.append(LineString(clean_arc(c)) if len(c) > 2 else a)
faces = list(polygonize(unary_union(new_arcs)))
# 2. faces back to parcels: the parcel that holds most of the face
tree = shapely.STRtree(orig)
owner = {}
for f in faces:
    cand = tree.query(f, predicate="intersects")
    best, ba = None, 0.0
    for k in cand:
        ar = f.intersection(orig[k]).area
        if ar > ba:
            best, ba = k, ar
    if best is not None and ba > 0.5 * f.area:
        owner.setdefault(best, []).append(f)
new = []
for k in range(len(orig)):
    parts = owner.get(k, [])
    g_ = unary_union(parts) if parts else orig[k]
    if g_.geom_type != "Polygon":
        g_ = max(getattr(g_, "geoms", [g_]), key=lambda p: p.area)
    new.append(g_)
# 3. back to EPSG:4326, untouched vertices to their exact original coordinates
ll = gpd.GeoSeries(new, crs=32644).to_crs(4326)
pts = np.vstack([shapely.get_coordinates(g_) for g_ in A.geometry])
kd = cKDTree(pts)
def back(gm):
    def f(c):
        d, k = kd.query(c); c2 = c.copy(); near = d <= 1e-8; c2[near] = pts[k[near]]; return c2
    return shapely.transform(gm, f)
B = A.copy(); B["geometry"] = [back(g_) for g_ in ll]
# 3b. small holes (seam wedges and hairlines, < FILL_M2) join the bordering parcel with the longest border,
#     a red (automation) parcel before one of his; the two 2 m strips along the red 70A stay open (his decision)
FILL_M2 = 15.0
K2_ = 110540.0 * 108600.0
geoms = list(B.geometry)
for _round in range(3):
    Ua = unary_union(geoms); filled = 0
    for p in getattr(Ua, "geoms", [Ua]):
        for r in p.interiors:
            h = Polygon(r).difference(Ua)
            if h.is_empty or h.area * K2_ > FILL_M2:
                continue
            cand = [(B.layer.iat[i].endswith("_auto"), h.boundary.intersection(geoms[i].buffer(1e-9)).length, i)
                    for i in range(len(geoms)) if geoms[i].distance(h) < 1e-9]
            if not cand:
                continue
            _red, _len, i = max(cand)
            u = unary_union([geoms[i], h])
            if u.geom_type == "Polygon":
                geoms[i] = u; filled += 1
    if not filled:
        break
B["geometry"] = geoms
new = list(gpd.GeoSeries(geoms, crs=4326).to_crs(32644))
# 4. report
rep = []
for i in range(len(A)):
    o, nw = orig[i], new[i]
    rep.append({"layer": A.layer.iat[i], "pos": A.pos.iat[i], "plot_no": A.plot_no.iat[i],
                "vertices_before": len(o.exterior.coords) - 1, "vertices_after": len(nw.exterior.coords) - 1,
                "area_change_m2": round(nw.area - o.area, 2), "area_change_%": round(100 * (nw.area / o.area - 1), 3),
                "max_move_m": round(shapely.hausdorff_distance(o, nw), 2), "valid": nw.is_valid,
                "changed": not A.geometry.iat[i].equals_exact(B.geometry.iat[i], 0)})
R = pd.DataFrame(rep)
pd.set_option("display.width", 200)
print(R[R.changed].to_string(index=False))
print("faces %d, parcels %d, parcels without a face %d" % (len(faces), len(orig), sum(1 for k in range(len(orig)) if k not in owner)))
K2 = 110540.0 * 108600.0
g = B.geometry.values; t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); k = a_ < b_
ov = [(B.survey.iat[i], B.survey.iat[j], round(sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") * K2, 4))
      for i, j in zip(a_[k], b_[k]) if B.survey.iat[i] != B.survey.iat[j] for x in [shapely.intersection(g[i], g[j])]]
ov = [o for o in ov if o[2] > 0]
Ua = unary_union(g); holes = [round(Polygon(r).difference(Ua).area * K2, 2) for p in getattr(Ua, "geoms", [Ua]) for r in p.interiors]
print("after: overlaps", len(ov), ov[:8], "| holes", holes)
B.to_file(OUT, layer="plan", driver="GPKG"); R.to_csv(OUT.replace(".gpkg", "_changes.csv"), index=False)
