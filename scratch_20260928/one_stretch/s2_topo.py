"""One Stretch, step 2: seams of the merged plots.
Akash, 2026-09-28: no new boundaries, no shape changed to fill a gap; what cannot be done that way is left for hand work.
So only HAIRLINES are closed - an overlap or a gap whose widest point is under 5 cm:
  - hairline overlap: the plot that gives way takes the other plot's existing line (largest vertex move under 6 cm);
  - hairline gap: cut across at the corners of the plots beside it, each piece joins the plot it lies along;
  - plots of Manual_Georeference are never changed (they are clean and approved): their neighbours close onto them.
Everything wider is left exactly as it is and listed.
-> out/merged_clean.gpkg (EPSG:32644), out/topo_log.csv, out/left_for_hand.gpkg"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union, split
HAIR = 0.05; MOVE = 0.06; VTOL = 0.005
G = gpd.read_file("out/merged_raw.gpkg"); n = len(G); geo = [g.buffer(0) if not g.is_valid else g for g in G.geometry]; orig = list(geo)
frozen = (G.src == "Manual_Georeference").values; name = lambda i: "%s %s/%s" % (G.village[i], G.sheet[i], G.first_read[i] or "-")
polys = lambda g: [q for q in shapely.get_parts(shapely.make_valid(g)) if q.geom_type == "Polygon" and q.area > 0]
one = lambda g: max(polys(g), key=lambda q: q.area)
wmax = lambda g: 2 * shapely.maximum_inscribed_circle(g).length if g.area > 0 else 0.0
log = []


def despike(g):
    c = [tuple(x) for x in np.array(g.exterior.coords)[:-1]]; ch = True
    while ch and len(c) > 3:
        ch = False; m = len(c)
        for i in range(m):
            a, b, d = np.array(c[i - 1]), np.array(c[i]), np.array(c[(i + 1) % m]); ab, bd = np.hypot(*(b - a)), np.hypot(*(d - b))
            if ab < 1e-7 or (bd > 0 and (a - b) @ (d - b) / ab / bd > 0.999999): del c[i]; ch = True; break
    return shapely.Polygon(c, [list(r.coords) for r in g.interiors])


def accept(i, g, what):
    """a changed shape is taken only if it is one valid polygon, keeps its holes, and no vertex moved beyond the limit"""
    if g is None or g.is_empty: return False
    ps = polys(g)
    if not ps: return False
    g = max(ps, key=lambda q: q.area)
    if sum(q.area for q in ps) - g.area > 0.01 and len(ps) > 1: return False                      # it would fall in two
    if len(g.interiors) != len(geo[i].interiors): return False
    far = max(geo[i].boundary.distance(shapely.Point(p)) for p in shapely.get_coordinates(g)); back = max(g.boundary.distance(shapely.Point(p)) for p in shapely.get_coordinates(geo[i]))
    if max(far, back) > MOVE: return False
    geo[i] = g; return True


tree = shapely.STRtree(geo); nb = {i: [j for j in tree.query(geo[i].buffer(0.1)) if j != i] for i in range(n)}
# A. vertices a few millimetres apart become one vertex: the plot takes its neighbours' vertices (Manual plots are the reference, never the mover)
moved_a = 0
for i in sorted(range(n), key=lambda i: (not frozen[i], -geo[i].area)):
    if frozen[i] or not nb[i]: continue
    ref = shapely.MultiPoint(np.vstack([shapely.get_coordinates(geo[j]) for j in nb[i]])); g = shapely.snap(geo[i], ref, VTOL)
    if g.is_valid and g.geom_type == "Polygon" and abs(g.area - geo[i].area) < 0.5 and g.hausdorff_distance(geo[i]) <= VTOL * 1.5 and not g.equals_exact(geo[i], 0): geo[i] = g; moved_a += 1
print("A. plots whose vertices were matched to their neighbours (under 5 mm):", moved_a)
# B. hairline overlaps
tree = shapely.STRtree(geo); nB = 0; areaB = 0.0; skipB = []
for i in range(n):
    for j in tree.query(geo[i]):
        if j <= i: continue
        x = geo[i].intersection(geo[j])
        if x.is_empty or x.area <= 0: continue
        hair = [q for q in polys(x) if wmax(q) <= HAIR]
        if not hair: continue
        if frozen[i] and frozen[j]: skipB.append((i, j, sum(q.area for q in hair), "both plots are Manual_Georeference")); continue
        give, keep = (j, i) if frozen[i] else (i, j) if frozen[j] else ((i, j) if geo[i].area >= geo[j].area else (j, i))
        h = unary_union(hair); ok = False
        for g_, k_ in ((give, keep), (keep, give)):
            if frozen[g_]: continue
            t = geo[g_].difference(h.buffer(1e-7, join_style=2)) if False else geo[g_].difference(h)
            if accept(g_, t, "overlap"): ok = True; log.append({"what": "hairline overlap closed", "plot": name(g_), "onto": name(k_), "sqm": round(h.area, 4), "widest_mm": round(1000 * max(wmax(q) for q in hair), 1)}); nB += 1; areaB += h.area; break
        if not ok: skipB.append((i, j, h.area, "the trimmed plot would break or move beyond 6 cm"))
print("B. hairline overlaps closed: %d (%.2f m2) | not closed: %d" % (nB, areaB, len(skipB)))
# C. hairline gaps
def attempt(h, rev):
    near = [k for k in tree.query(h.buffer(0.01)) if geo[k].distance(h) < 1e-6]; movers = [k for k in near if not frozen[k]]
    trial = {}; done = []; rest = [h]
    lies = lambda x, k: max(geo[k].distance(shapely.Point(p)) for p in shapely.get_coordinates(x)) < MOVE
    cur = lambda k: trial.get(k, geo[k]); share = lambda k, piece: cur(k).boundary.intersection(piece.boundary.buffer(1e-5)).length
    order = sorted(movers, key=lambda k: -geo[k].boundary.intersection(h.boundary.buffer(1e-5)).length); order = order[::-1] if rev else order
    for k in order:
        new_rest = []
        for piece in rest:
            seg = geo[k].boundary.intersection(piece.boundary.buffer(1e-6)); cuts = []
            if seg.length >= 1e-4:
                seg = shapely.line_merge(seg) if seg.geom_type == "MultiLineString" else seg
                for ln in getattr(seg, "geoms", [seg]):
                    cc = np.array(ln.coords)
                    if len(cc) < 2: continue
                    for p, o in ((cc[0], cc[1]), (cc[-1], cc[-2])):
                        if np.hypot(*(o - p)) < 1e-9: continue
                        t = (o - p) / np.hypot(*(o - p)); nrm = np.array([-t[1], t[0]]); cuts.append(shapely.LineString([p - nrm * (HAIR * 2), p + nrm * (HAIR * 2)]))
            sub = [piece]
            for c in cuts:
                nxt = []
                for x in sub:
                    try: nxt += [y for y in split(x, c).geoms if y.area > 0]
                    except Exception: nxt.append(x)
                sub = nxt
            for x in sub:
                if lies(x, k) and share(k, x) > 0.35 * x.length / 2 and share(k, x) > 1e-4:
                    g = shapely.make_valid(unary_union([cur(k), x.buffer(1e-7, join_style=2)])); ps = polys(g)
                    if len(ps) == 1 and len(ps[0].interiors) == len(geo[k].interiors): trial[k] = despike(ps[0]); done.append((k, x.area)); continue
                new_rest.append(x)
        rest = new_rest
        if not rest: break
    for x in list(rest):
        cand = sorted(((share(k, x), k) for k in movers if lies(x, k)), reverse=True)
        if cand and cand[0][0] > 1e-4:
            k = cand[0][1]; g = shapely.make_valid(unary_union([cur(k), x.buffer(1e-7, join_style=2)])); ps = polys(g)
            if len(ps) == 1 and len(ps[0].interiors) == len(geo[k].interiors): trial[k] = despike(ps[0]); done.append((k, x.area)); rest.remove(x)
    left = sum(x.area for x in rest); ok = left < 0.0005 or left < 0.02 * h.area
    for k, g in trial.items():
        if not g.is_valid or geo[k].symmetric_difference(g).area > h.area + 0.002 or geo[k].hausdorff_distance(g) > MOVE: ok = False
    return ok, trial, done, left, near


nC = 0; areaC = 0.0; skipC = []
for rnd in range(3):
    tree = shapely.STRtree(geo); U = unary_union(geo); holes = [shapely.Polygon(r) for q in polys(U) for r in q.interiors]; holes = [h for h in holes if h.area > 0 and wmax(h) <= HAIR]
    if not holes: break
    closed = 0
    for h in sorted(holes, key=lambda h: -h.area):
        for rev in (False, True):
            ok, trial, done, left, near = attempt(h, rev)
            if ok and trial: break
        if ok and trial:
            for k, g in trial.items(): geo[k] = g
            nC += 1; areaC += h.area; closed += 1; log.append({"what": "hairline gap closed", "plot": ", ".join(sorted({name(k) for k, _a in done})), "onto": "", "sqm": round(h.area, 4), "widest_mm": round(1000 * wmax(h), 1)})
        elif rnd == 2: skipC.append((h, "its pieces would leave the plots' own lines or only Manual plots lie beside it"))
    # a fill may touch a neighbour by a hair: the filled plot gives that hair back
    tree = shapely.STRtree(geo)
    for i in range(n):
        if frozen[i]: continue
        for j in tree.query(geo[i]):
            if j != i:
                x = geo[i].intersection(geo[j])
                if 0 < x.area < 1e-4 and (frozen[j] or j > i): accept(i, geo[i].difference(geo[j]), "micro")
    if not closed: break
print("C. hairline gaps closed: %d (%.2f m2) | not closed: %d" % (nC, areaC, len(skipC)))
# what stays for hand work, and the proof that nothing else changed
tree = shapely.STRtree(geo); left = []
for i in range(n):
    for j in tree.query(geo[i]):
        if j > i:
            x = geo[i].intersection(geo[j])
            for q in polys(x) if x.area > 0 else []:
                if q.area > 1e-6: left.append({"kind": "overlap", "between": name(i) + " | " + name(j), "layers": G.src[i] + " | " + G.src[j], "sqm": round(q.area, 3), "widest_m": round(wmax(q), 3), "geometry": q})
U = unary_union(geo)
for q in polys(U):
    for r in q.interiors:
        h = shapely.Polygon(r)
        if h.area > 1e-6:
            near = sorted({name(k) for k in tree.query(h.buffer(0.01)) if geo[k].distance(h) < 1e-6}); left.append({"kind": "gap", "between": " | ".join(near), "layers": " | ".join(sorted({G.src[k] for k in tree.query(h.buffer(0.01)) if geo[k].distance(h) < 1e-6})), "sqm": round(h.area, 3), "widest_m": round(wmax(h), 3), "geometry": h})
L = gpd.GeoDataFrame(left, geometry="geometry", crs=32644); L.to_file("out/left_for_hand.gpkg", driver="GPKG")
ch = [i for i in range(n) if not geo[i].equals_exact(orig[i], 0)]; mv = [orig[i].hausdorff_distance(geo[i]) for i in ch]; da = [abs(geo[i].area - orig[i].area) for i in ch]
print("plots changed: %d of %d | Manual plots changed: %d | largest vertex move %.1f mm | largest area change %.3f m2 | total area %.2f -> %.2f m2" % (len(ch), n, sum(1 for i in ch if frozen[i]), 1000 * max(mv or [0]), max(da or [0]), sum(g.area for g in orig), sum(g.area for g in geo)))
assert not any(frozen[i] for i in ch) and max(mv or [0]) <= MOVE + 1e-9 and all(g.is_valid and g.geom_type == "Polygon" for g in geo)
print("left for hand work: overlaps %d (%.1f m2), gaps %d (%.1f m2)" % ((L.kind == "overlap").sum(), L[L.kind == "overlap"].sqm.sum(), (L.kind == "gap").sum(), L[L.kind == "gap"].sqm.sum()))
print("   of these under 5 cm at the widest (could not be closed without breaking a rule):", int((L.widest_m <= HAIR).sum()))
G2 = G.copy(); G2["geometry"] = geo; G2["moved_mm"] = [round(1000 * orig[i].hausdorff_distance(geo[i]), 2) for i in range(n)]; G2.to_file("out/merged_clean.gpkg", driver="GPKG"); pd.DataFrame(log).to_csv("out/topo_log.csv", index=False)
