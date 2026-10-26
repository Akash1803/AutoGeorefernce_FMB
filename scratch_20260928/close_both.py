"""Village 35_04_074 seam clean-up (Akash 2026-09-28: keep 47B, close the gaps without changing shapes).
1. overlaps: the loser plot is trimmed by the keeper (582 / 581 / 47B always keep; otherwise the larger plot gives way)
2. enclosed hairline holes: each is cut across at the corners of the receiving plots and every piece joins the plot it
   lies along, so a plot only moves its own edge by the width of the hairline.
Dry run by default: writes out/close_plan.json {file: {fid: wkt}} and prints what changes."""
import warnings, sys, json, os; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union, split
from autogeoref import visible, paths, fabric_check as fc
src = dict(visible.shown("35_04_074")); src.update({"K" + k: v for k, v in visible.shown("35_04_077").items()})
KEEP = ["582", "581", "47B"]; OV_MIN = 0.001; HOLE_MIN = float(os.environ.get("HOLE_MIN", "0.001"))
parts = []
for s, (f, l) in src.items():
    g = gpd.read_file(f, layer=l, fid_as_index=True); assert g.crs.to_epsg() == 32644, (s, g.crs)
    g = g[g.geometry.notna() & ~g.geometry.is_empty].copy(); g["fid"] = g.index; g["s"] = s; g["file"] = str(f); g["layer"] = l
    parts.append(g[["s", "plot_no", "fid", "file", "layer", "geometry"]])
P = gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=32644)
P["multi"] = [x.geom_type.startswith("Multi") for x in P.geometry]
P["geometry"] = [max(fc._polys(x.buffer(0)), key=lambda q: q.area) for x in P.geometry]; P["orig"] = P.geometry.copy()
name = lambda i: f"{P.s[i]}/{P.plot_no[i]}"; log = []; changed = set()
one = lambda g: max(fc._polys(shapely.make_valid(g)), key=lambda q: q.area)
# ---- 1. overlaps
sidx = P.sindex
for i in range(len(P)):
    for j in sidx.query(P.geometry[i]):
        if j <= i: continue
        a = P.geometry[i].intersection(P.geometry[j]).area
        if a < OV_MIN: continue
        ki, kj = P.s[i] in KEEP, P.s[j] in KEEP
        if ki and kj: rank = lambda k: KEEP.index(P.s[k]); lose, keep = (j, i) if rank(i) < rank(j) else (i, j)
        elif ki: lose, keep = j, i
        elif kj: lose, keep = i, j
        else: lose, keep = (i, j) if P.geometry[i].area >= P.geometry[j].area else (j, i)
        t = shapely.make_valid(P.geometry[lose].difference(P.geometry[keep])); ps = sorted(fc._polys(t), key=lambda q: -q.area)
        log.append(f"overlap {a:.3f} m2: {name(lose)} gives way to {name(keep)}" + (f" (dropped pieces {[round(q.area, 3) for q in ps[1:]]})" if len(ps) > 1 else ""))
        P.at[lose, "geometry"] = ps[0]; changed.add(lose)
# ---- 2. holes
def despike(g):
    c = [tuple(x) for x in np.array(g.exterior.coords)[:-1]]; ch = True
    while ch and len(c) > 3:
        ch = False; n = len(c)
        for i in range(n):
            a, b, d = np.array(c[i - 1]), np.array(c[i]), np.array(c[(i + 1) % n]); ab, bd = np.hypot(*(b - a)), np.hypot(*(d - b))
            if ab < 1e-7 or (bd > 0 and (a - b) @ (d - b) / ab / bd > 0.999999): del c[i]; ch = True; break
    return shapely.Polygon(c, [list(r.coords) for r in g.interiors])
geo = unary_union(list(P.geometry)); holes = [shapely.Polygon(r) for q in fc._polys(geo) for r in q.interiors if HOLE_MIN <= shapely.Polygon(r).area < 50]
filled = unary_union([shapely.Polygon(q.exterior) for q in fc._polys(geo)])
for q in fc._polys(filled.buffer(0.06, join_style=2).buffer(-0.06, join_style=2).difference(filled)):       # hairlines that reach the outer edge
    if q.area < HOLE_MIN or 2 * q.area / q.length > 0.1: continue
    nb = [k for k in P.sindex.query(q.buffer(0.01)) if P.geometry[k].boundary.intersection(q.boundary.buffer(1e-5)).length > 0.5]
    if len({P.s[k] for k in nb}) >= 2: holes.append(q); print(f"open hairline {q.area:.3f} m2 at ({q.centroid.x:.0f},{q.centroid.y:.0f}), mean width {2000 * q.area / q.length:.1f} mm, between {sorted({P.s[k] for k in nb})}")
holes = sorted(holes, key=lambda h: -h.area)
sidx = P.sindex; skipped = []
def attempt(h, rev):
    near = [k for k in sidx.query(h.buffer(0.01)) if P.geometry[k].distance(h) < 1e-6]
    tag = f"hole {h.area:.3f} m2 at ({h.centroid.x:.0f},{h.centroid.y:.0f})"; wmax = 2 * shapely.maximum_inscribed_circle(h).length if h.area > 0 else 0
    rest = [h]; trial = {}; done = []
    lies = lambda x, k: max(P.geometry[k].distance(shapely.Point(p)) for p in shapely.get_coordinates(x)) < 0.15
    share = lambda k, piece: P.geometry[k].boundary.intersection(piece.boundary.buffer(1e-5)).length if k not in trial else trial[k].boundary.intersection(piece.boundary.buffer(1e-5)).length
    surveys_ = sorted({P.s[k] for k in near if P.s[k] not in KEEP[:2]}, key=lambda s: -sum(P.geometry[k].boundary.intersection(h.boundary.buffer(1e-5)).length for k in near if P.s[k] == s))
    surveys = surveys_[::-1] if rev else surveys_
    for s in surveys:
        mine = [k for k in near if P.s[k] == s]; new_rest = []
        for piece in rest:
            cuts = []; seen = set()
            for k in mine:
                seg = P.geometry[k].boundary.intersection(piece.boundary.buffer(1e-6))
                if seg.length < 1e-4: continue
                seg = shapely.line_merge(seg) if seg.geom_type == "MultiLineString" else seg
                for ln in getattr(seg, "geoms", [seg]):
                    cc = np.array(ln.coords)
                    if len(cc) < 2: continue
                    for p, o in ((cc[0], cc[1]), (cc[-1], cc[-2])):
                        key = (round(p[0], 4), round(p[1], 4))
                        if key in seen or np.hypot(*(o - p)) < 1e-9: continue
                        seen.add(key); t = (o - p) / np.hypot(*(o - p)); nrm = np.array([-t[1], t[0]])
                        cuts.append(shapely.LineString([p - nrm * 0.2, p + nrm * 0.2]))
            sub = [piece]
            for c in cuts:
                nxt = []
                for x in sub:
                    try: nxt += [y for y in split(x, c).geoms if y.area > 0]
                    except Exception: nxt.append(x)
                sub = nxt
            for x in sub:
                cand = sorted(((share(k, x), k) for k in mine if lies(x, k)), reverse=True)
                if cand and cand[0][0] > 0.35 * x.length / 2 and cand[0][0] > 1e-4:
                    k = cand[0][1]; g = shapely.make_valid(unary_union([trial.get(k, P.geometry[k]), x.buffer(1e-7, join_style=2)]))
                    if len(fc._polys(g)) == 1 and not fc._polys(g)[0].interiors: trial[k] = despike(fc._polys(g)[0]); done.append((k, x.area)); continue
                new_rest.append(x)
        rest = new_rest
        if not rest: break
    for x in list(rest):                                                                  # compact leftovers: to the plot they lie along most
        cand = sorted(((share(k, x), k) for k in near if P.s[k] not in KEEP[:2] and lies(x, k)), reverse=True)
        if cand and cand[0][0] > 1e-4:
            k = cand[0][1]; g = shapely.make_valid(unary_union([trial.get(k, P.geometry[k]), x.buffer(1e-7, join_style=2)]))
            if len(fc._polys(g)) == 1 and not fc._polys(g)[0].interiors: trial[k] = despike(fc._polys(g)[0]); done.append((k, x.area)); rest.remove(x)
    left = sum(x.area for x in rest); ok = left < 0.0005 or left < 0.02 * h.area
    for k, g in trial.items():
        o = P.geometry[k]; hd = o.hausdorff_distance(g)
        if not g.is_valid or o.symmetric_difference(g).area > h.area + 0.002 or (hd > max(0.12, 3 * wmax) and o.symmetric_difference(o.buffer(0)).area > 1e-9): ok = False; tag += f" [{name(k)}: area {g.area - o.area:+.3f}, moved {hd:.3f} m]"
    return ok, trial, done, tag, wmax, left, near
for h in holes:
    for rev in (False, True):
        ok, trial, done, tag, wmax, left, near = attempt(h, rev)
        if ok: break
    if ok:
        for k, g in trial.items(): P.at[k, "geometry"] = g; changed.add(k)
        log.append(f"{tag}, width up to {1000 * wmax:.1f} mm -> " + ", ".join(f"{name(k)} +{a:.4f}" for k, a in done))
    else: skipped.append(f"{tag}, width up to {1000 * wmax:.1f} mm, between {sorted({P.s[k] for k in near})} (left {left:.4f})")
# micro overlaps the fills may have made against kept plots
for k in list(changed):
    for j in sidx.query(P.geometry[k]):
        if j != k and P.geometry[k].intersection(P.geometry[j]).area > 1e-7 and (j not in changed or P.s[j] in KEEP):
            P.at[k, "geometry"] = one(P.geometry[k].difference(P.geometry[j]))
print("\n".join(log)); print("\nNOT closed:", *skipped, sep="\n   ")
plan = {}
print("\nplots changed:", len(changed))
for k in sorted(changed, key=lambda k: (P.s[k], str(P.plot_no[k]))):
    o, g = P.orig[k], P.geometry[k]; assert g.is_valid and g.geom_type == "Polygon"
    far = max(o.distance(shapely.Point(p)) for p in shapely.get_coordinates(g)); assert far < 0.2 or o.symmetric_difference(g).area > 0.5, (name(k), far)
    print(f"   {name(k):12s} {o.area:10.3f} -> {g.area:10.3f} m2 ({g.area - o.area:+.3f}) | land changed {o.symmetric_difference(g).area:.3f} m2 | largest move {o.hausdorff_distance(g):.3f} m | vertices {len(o.exterior.coords) - 1} -> {len(g.exterior.coords) - 1}")
    plan.setdefault(P.file[k] + "|" + P.layer[k], {})[str(int(P.fid[k]))] = g.wkt
json.dump(plan, open("out/close_plan2.json", "w")); print("files:", len(plan))
P2 = P.copy(); sidx = P2.sindex; ov = []
for i in range(len(P2)):
    for j in sidx.query(P2.geometry[i]):
        if j > i and P2.geometry[i].intersection(P2.geometry[j]).area > 0.001: ov.append((name(i), name(j), round(P2.geometry[i].intersection(P2.geometry[j]).area, 4)))
g2 = unary_union(list(P2.geometry)); hs = [shapely.Polygon(r).area for q in fc._polys(g2) for r in q.interiors]
print("AFTER: overlaps > 0.001 m2:", ov, "| holes > 0.001 m2:", [round(a, 4) for a in hs if a > 0.001], "| all holes:", len(hs), "total", round(sum(hs), 4))
