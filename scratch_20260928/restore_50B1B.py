"""Put 50B/1B back into Manual_Georeference as the exact space left between the plots Akash fixed, and close the two hairlines
his edit left along 582. Nothing is moved: the new plot takes its neighbours' own vertices; the hairlines are closed by giving
40A and 40A/2A1 the vertex 582 already has there (0.86 mm and 0.15 mm) and 582 one node on its own straight edge.
DRY=1 only reports."""
import warnings, os, sys, json; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import geopandas as gpd, shapely, numpy as np
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union
from autogeoref import qgis_bridge as qb
LIVE = r"D:/Projects/Tambaram_Chengalpattu/Manual_Georeferenced/Manual_Georeference.geojson"
B = gpd.read_file(LIVE); M = B.to_crs(32644); old = gpd.read_file("out/rebuild/Manual_Georeference.geojson"); old = old[old.kide == "50B/1B"].iloc[0]
assert len(B) == 661 and not (B.kide == "50B/1B").any()
U = unary_union(list(B.geometry)); holes = sorted([Polygon(r) for p in (U.geoms if U.geom_type == "MultiPolygon" else [U]) for r in p.interiors], key=lambda h: -h.area)
to_m = lambda g: gpd.GeoSeries([g], crs=B.crs).to_crs(32644).iloc[0]
big = holes[0]; bm = to_m(big); assert 430 < bm.area < 438 and bm.intersection(to_m(old.geometry)).area / bm.area > 0.99, bm.area
own = {c for g in B.geometry for c in g.exterior.coords}; assert all(c in own for c in big.exterior.coords), "a corner of the space is not a corner of a neighbour"
new = shapely.geometry.polygon.orient(big, 1.0)
# hairlines: every corner of a hairline that one neighbour has and the other only passes -> give it to the one that passes
ins = []   # (kide, area_sqm of the plot, lon, lat)
for h in holes[1:]:
    hm = to_m(h); assert hm.area < 0.02, hm.area
    for c in list(h.exterior.coords)[:-1]:
        for j in B.index[M.geometry.distance(to_m(Point(c))) < 0.002]:
            if c not in set(B.geometry[j].exterior.coords): ins.append((B.kide[j], float(B.area_sqm[j]), j, c))
print("new 50B/1B: %.1f m2, %d corners" % (bm.area, len(new.exterior.coords) - 1)); [print("   node for", k, a, "at", c, "| %.2f mm from its edge" % (M.geometry[j].exterior.distance(to_m(Point(c))) * 1000)) for k, a, j, c in ins]
assert sorted(k for k, *_ in ins) == ["40A", "40A/2A1", "582"], ins
def with_node(g, c):
    r = list(g.exterior.coords); m = to_m(g); p = to_m(Point(c)); best = None
    for i in range(len(r) - 1):
        d = shapely.LineString([r[i], r[i + 1]]).distance(Point(c))
        if best is None or d < best[0]: best = (d, i)
    r.insert(best[1] + 1, c); q = Polygon(r, [list(x.coords) for x in g.interiors]); assert q.is_valid and abs(to_m(q).area - m.area) < 0.02; return q
geo = {}
for k, a, j, c in ins: geo[j] = with_node(geo.get(j, B.geometry[j]), c)
# proof before anything is written
T = B.copy(); 
for j, g in geo.items(): T.at[j, "geometry"] = g
allg = list(T.geometry) + [new]; U2 = unary_union(allg); P2 = list(U2.geoms) if U2.geom_type == "MultiPolygon" else [U2]
h2 = [r for p in P2 for r in p.interiors]; tr = shapely.STRtree(allg); ov = sum(1 for i, g in enumerate(allg) for j in tr.query(g) if j > i and g.intersection(allg[j]).area > 1e-16)
print("after: pieces", len(P2), "| holes", len(h2), "| overlaps", ov); assert len(P2) == 1 and not h2 and ov == 0
# areas true to the shapes he left
changed = {j: round(to_m(T.geometry[j]).area, 1) for j in B.index if abs(to_m(B.geometry[j]).area - float(B.area_sqm[j])) > 0.06}
print("area fields to refresh:", {B.kide[j]: (float(B.area_sqm[j]), v) for j, v in changed.items()})
att = {c: old[c] for c in B.columns if c != "geometry"}; att["area_sqm"] = round(bm.area, 1); att["area_acre"] = round(bm.area / 4046.8564224, 4); print(att)
pay = {"new_wkb": new.wkb_hex, "att": att, "geo": [{"kide": B.kide[j], "area_sqm": float(B.area_sqm[j]), "wkb": g.wkb_hex} for j, g in geo.items()], "areas": [{"kide": B.kide[j], "old": float(B.area_sqm[j]), "sqm": v, "acre": round(v / 4046.8564224, 4)} for j, v in changed.items()]}
json.dump(pay, open("out/restore_50B1B.json", "w"), indent=1)
if os.environ.get("DRY"): sys.exit()
code = r'''
import json
from qgis.core import QgsProject, QgsGeometry, QgsFeature
p = json.load(open(PAY)); l = QgsProject.instance().mapLayersByName("Manual_Georeference")[0]
assert not l.isEditable() and not l.isModified(), "the layer is being edited"
assert l.featureCount() == 661 and not [f for f in l.getFeatures() if f["kide"] == "50B/1B"]
idx = {f.name(): i for i, f in enumerate(l.fields())}; feats = list(l.getFeatures())
def one(k, a):
    m = [f for f in feats if str(f["kide"]) == k and abs(float(f["area_sqm"]) - a) < 0.05]; assert len(m) == 1, (k, a, len(m)); return m[0]
def geom(h):
    g = QgsGeometry(); g.fromWkb(bytes.fromhex(h)); assert not g.isNull() and g.isGeosValid(); return g
G = [(one(g["kide"], g["area_sqm"]).id(), geom(g["wkb"])) for g in p["geo"]]; A = [(one(a["kide"], a["old"]).id(), a) for a in p["areas"]]; N = geom(p["new_wkb"])
l.startEditing()
try:
    for fid, g in G: assert l.changeGeometry(fid, g)
    for fid, a in A: assert l.changeAttributeValue(fid, idx["area_sqm"], a["sqm"]) and l.changeAttributeValue(fid, idx["area_acre"], a["acre"])
    f = QgsFeature(l.fields()); f.setGeometry(N)
    for k, v in p["att"].items(): f.setAttribute(idx[k], v)
    assert l.addFeature(f); ok = l.commitChanges()
    if not ok: l.rollBack()
except Exception:
    l.rollBack(); raise
l.triggerRepaint(); print(json.dumps({"committed": ok, "errors": l.commitErrors()[:3], "features": l.featureCount()}))
'''.replace("PAY", repr(os.path.abspath("out/restore_50B1B.json").replace("\\", "/")))
print(qb.json_result(code, timeout=200))
