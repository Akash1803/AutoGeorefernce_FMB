"""Manual_Georeference.geojson for Thirukatchur + Kizhikaranai from the layers QGIS shows, as a clean coverage:
every seam is one line with the same vertices on both sides, so the QGIS Topology Checker (gaps / duplicates /
invalid / multi-part) finds nothing. Seam work here is sub-millimetre: vertices within 3 mm are made common, leftover
micro overlaps trimmed, micro holes given to the plot they lie along, neighbour vertices inserted on shared edges."""
import warnings, sys, os, re, shutil, datetime; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
from autogeoref import visible, fabric_check as fc
OUT = r"D:\Projects\Tambaram_Chengalpattu\Manual_Georeferenced\Manual_Georeference.geojson"; TOL = 0.003
BUF = unary_union(list(gpd.read_file(r"D:\Projects\Tambaram_Chengalpattu\Tambaram_Chengalpattu_Railway.gpkg", layer="rail_buffer_30m").to_crs(32644).geometry))


def status(src, single, pn):
    src = "" if src is None else str(src)
    if "a sliver left by the placement" in src: return "sliver"
    if single: return "whole survey - the FMB sheet has no subdivision"
    if not pn: return "no number printed on the FMB sheet" + (" - the plot lies outside the survey's recorded area" if "outside the recorded area" in src else "")
    if "no number printed" in src and "missing" in src: return "not printed on the FMB sheet - taken as the only number missing on it"
    return "printed on the FMB sheet"
rows = []; left_out = []
key = lambda s: (int(re.match(r"\d+", s).group(0)), s)
for V, name in [("35_04_074", "Thirukatchur"), ("35_04_077", "Kizhikaranai")]:
    sh = visible.shown(V)
    for s in sorted(sh, key=key):
        g = gpd.read_file(sh[s][0], layer=sh[s][1]); g = g.to_crs(32644) if g.crs.to_epsg() != 32644 else g
        single = int((g.geometry.notna() & ~g.geometry.is_empty).sum()) == 1
        for _, r in g.iterrows():
            if r.geometry is None or r.geometry.is_empty: left_out.append((name, s, "no geometry")); continue
            geom = max(fc._polys(r.geometry.buffer(0)), key=lambda q: q.area); pn = r.get("plot_no"); pn = "" if pn is None or str(pn) in ("nan", "None", "NULL", "0", "") else str(pn).strip()
            if single or pn == s: pn = ""
            st = status(r.get("label_source"), single, pn)
            rows.append({"survey_no": s, "subdiv_no": pn, "village": name, "village_code": V, "subdiv_status": st, "geometry": geom})
# slivers of the placement are not plots: each joins the neighbour it shares most of its edge with
keep = []
for i, r in enumerate(rows):
    if r["geometry"].area < 0.5: left_out.append((r["village"], r["survey_no"], "sliver %.2f m2" % r["geometry"].area)); continue      # closed by the seam clean-up below
    if r["subdiv_status"] == "sliver":
        nb = sorted(((q["geometry"].boundary.intersection(r["geometry"].buffer(0.01)).length, j) for j, q in enumerate(rows) if j != i and q["subdiv_status"] != "sliver" and q["geometry"].area >= 0.5 and q["geometry"].distance(r["geometry"]) < 0.01), reverse=True)
        done = False
        for _, j in nb[:3]:
            u = shapely.make_valid(unary_union([rows[j]["geometry"], shapely.snap(r["geometry"], rows[j]["geometry"], 0.003)])); ps = fc._polys(u)
            if len(ps) == 1 and not ps[0].interiors:
                rows[j]["geometry"] = ps[0]; left_out.append((r["village"], r["survey_no"], "sliver %.2f m2 joined to %s/%s" % (r["geometry"].area, rows[j]["survey_no"], rows[j]["subdiv_no"]))); done = True; break
        if done: continue
    keep.append(i)
rows = [rows[i] for i in keep]
G = gpd.GeoDataFrame(rows, geometry="geometry", crs=32644); n = len(G); geo = list(G.geometry); a0 = [g.area for g in geo]
one = lambda g: max(fc._polys(shapely.make_valid(g)), key=lambda q: q.area)


def despike(g):
    c = [tuple(x) for x in np.array(g.exterior.coords)[:-1]]; ch = True
    while ch and len(c) > 3:
        ch = False; m = len(c)
        for i in range(m):
            a, b, d = np.array(c[i - 1]), np.array(c[i]), np.array(c[(i + 1) % m]); ab, bd = np.hypot(*(b - a)), np.hypot(*(d - b))
            if ab < 1e-6 or (bd > 0 and (a - b) @ (d - b) / ab / bd > 0.999999): del c[i]; ch = True; break
    return shapely.Polygon(c, [list(r.coords) for r in g.interiors])


tree = shapely.STRtree(geo); nb = {i: [j for j in tree.query(geo[i].buffer(0.05)) if j != i] for i in range(n)}
# 1. common vertices: snap each plot to the vertices of its neighbours (already-snapped ones first)
for i in range(n):
    if not nb[i]: continue
    ref = shapely.MultiPoint(np.vstack([shapely.get_coordinates(geo[j]) for j in nb[i]]))
    g = shapely.snap(geo[i], ref, TOL); g = one(g) if not g.is_valid or g.geom_type != "Polygon" else g
    if abs(g.area - geo[i].area) < 0.5: geo[i] = despike(g)
# 2. micro overlaps: the larger plot gives way
for i in range(n):
    for j in nb[i]:
        if j > i and geo[i].intersection(geo[j]).area > 1e-9:
            k, o = (i, j) if geo[i].area >= geo[j].area else (j, i); geo[k] = one(geo[k].difference(geo[o]))
# 3-5. repeat until stable: holes to the plot they lie along, seams noded, coordinates on a 0.001 mm grid, overlaps trimmed
def node(ring, allv, vt):
    c = np.array(ring.coords); out = []
    for a, b in zip(c[:-1], c[1:]):
        out.append(tuple(a)); seg = shapely.LineString([a, b])
        if seg.length < 1e-9: continue
        mids = [allv[k] for k in vt.query(seg.buffer(1e-3)) if 5e-4 < np.hypot(*(allv[k] - a)) and 5e-4 < np.hypot(*(allv[k] - b)) and seg.distance(shapely.Point(allv[k])) < 6e-4]
        out += [tuple(v) for v in sorted(mids, key=lambda v: np.hypot(*(v - a)))]
    return out
for rnd in range(8):
    geo = [g if g.is_valid and g.geom_type == "Polygon" else one(g) for g in geo]
    tree = shapely.STRtree(geo); nov = 0
    for i in range(n):
        for j in tree.query(geo[i]):
            if j > i and geo[i].intersection(geo[j]).area > 0:
                k, o = (i, j) if geo[i].area >= geo[j].area else (j, i); geo[k] = one(geo[k].difference(geo[o])); nov += 1
    geo = [shapely.Polygon(g.exterior, [r for r in g.interiors if shapely.Polygon(r).area > 1.0]) for g in geo]      # a micro ring inside one plot is not a real hole
    u = unary_union(geo); holes = [shapely.Polygon(r) for q in fc._polys(u) for r in q.interiors]; tree = shapely.STRtree(geo)
    for h in holes:
        assert h.area < 0.3, ("unexpected gap", h.area, h.centroid.wkt)
        cand = sorted(((geo[k].boundary.intersection(h.boundary.buffer(1e-6)).length, k) for k in tree.query(h.buffer(0.01))), reverse=True)
        for _, k in cand:
            g = shapely.make_valid(unary_union([geo[k], h.buffer(1e-8, join_style=2)]))
            if len(fc._polys(g)) == 1 and len(fc._polys(g)[0].interiors) == len(geo[k].interiors): geo[k] = despike(fc._polys(g)[0]); break
    allv = np.unique(np.round(np.vstack([shapely.get_coordinates(g) for g in geo]), 3), axis=0); vt = shapely.STRtree(shapely.points(allv))
    for i in range(n):
        g = shapely.Polygon(node(geo[i].exterior, allv, vt), [node(r, allv, vt) for r in geo[i].interiors]); geo[i] = g if g.is_valid else geo[i]
    geo = [shapely.transform(g, lambda c: np.round(c, 3)) for g in geo]
    geo = [despike(g) if len(g.exterior.coords) > 4 else g for g in geo]
    bad = sum(1 for g in geo if not g.is_valid); print(f"round {rnd}: overlaps trimmed {nov}, holes filled {len(holes)}, invalid after rounding {bad}")
    if rnd > 0 and not holes and bad == 0 and nov <= 6: break
tree = shapely.STRtree(geo)
for i in range(n):                                   # last pass: the larger plot takes the exact edge of the smaller one, vertex for vertex
    for j in tree.query(geo[i]):
        if j > i and geo[i].intersection(geo[j]).area > 0:
            k, o = (i, j) if geo[i].area >= geo[j].area else (j, i)
            ref = shapely.MultiPoint(shapely.get_coordinates(geo[o])); g = shapely.snap(geo[k], ref, 0.0015)
            g = one(g.difference(geo[o])); g = shapely.transform(g, lambda c: np.round(c, 3)); geo[k] = despike(g) if g.is_valid else geo[k]
G["geometry"] = geo; assert all(g.is_valid and g.geom_type == "Polygon" for g in geo)
moved = [a.hausdorff_distance(b) for a, b in zip([r["geometry"] for r in rows], geo)]
G["kide"] = [f"{s}/{p}" if p else s for s, p in zip(G.survey_no, G.subdiv_no)]; G["area_sqm"] = [round(g.area, 1) for g in geo]; G["area_acre"] = [round(g.area / 4046.8564224, 4) for g in geo]
G["taluk"] = "Chengalpattu"; G["district"] = "Chengalpattu"; G["status"] = "QC Passed"; G["source_pdf"] = ["FMB_Sketches/%s/%s.pdf" % (v, s) for v, s in zip(G.village_code, G.survey_no)]
G["in_rail_buffer_30m"] = ["yes" if g.intersects(BUF) else "no" for g in geo]
G = G[["kide", "survey_no", "subdiv_no", "area_acre", "area_sqm", "village", "village_code", "taluk", "district", "status", "subdiv_status", "in_rail_buffer_30m", "source_pdf", "geometry"]]
W = G.to_crs(4326); W["geometry"] = [shapely.transform(g, lambda c: np.round(c, 10)) for g in W.geometry]
if os.path.exists(OUT):
    bk = r"D:\Projects\Tambaram_Chengalpattu\_logs\Manual_Georeference_before_%s.geojson" % datetime.datetime.now().strftime("%Y%m%d_%H%M%S"); shutil.copy2(OUT, bk); print("old file kept as", bk)
tmp = OUT + ".tmp.geojson"; W.to_file(tmp, driver="GeoJSON", COORDINATE_PRECISION=10)
open(OUT, "w", encoding="utf-8").write(open(tmp, encoding="utf-8").read().replace('"name": "Manual_Georeference.geojson.tmp"', '"name": "Manual_Georeference"')); os.remove(tmp)
B = gpd.read_file(OUT); u = unary_union(list(B.geometry)); rings = sum(len(q.interiors) for q in getattr(u, "geoms", [u]))
sidx = B.sindex; ov = [(B.kide[i], B.kide[j]) for i, g in enumerate(B.geometry) for j in sidx.query(g) if j > i and g.intersection(B.geometry[j]).area > 1e-16]
dup = [(B.kide[i], B.kide[j]) for i, g in enumerate(B.geometry) for j in sidx.query(g) if j > i and g.equals(B.geometry[j])]
print("plots", len(B), "| left out:", left_out); print(B.groupby("village").agg(surveys=("survey_no", "nunique"), plots=("kide", "size"), acres=("area_acre", "sum")).round(2).to_string())
print("CHECK on the written file (EPSG:4326): gaps (interior rings)", rings, "| fabric pieces", len(getattr(u, "geoms", [u])), "| overlaps", len(ov), ov[:5], "| duplicates", dup, "| invalid", int((~B.geometry.is_valid).sum()), "| multi-part", int((B.geom_type != "Polygon").sum()))
print("largest move of any plot edge: %.4f m | plots moved more than 1 mm: %d | largest area change %.3f m2" % (max(moved), sum(1 for m in moved if m > 0.001), max(abs(g.area - a) for g, a in zip(geo, a0))))
