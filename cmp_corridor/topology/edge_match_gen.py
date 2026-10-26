"""Close the visible gaps of CMP_my_Part by EDGE MATCHING (no new polygons): for each gap, ONE side (the receiving parcels of one
survey) is extended onto the opposite boundary; where two receiving plots meet at the gap, their dividing line is carried straight
across (no Voronoi arcs, no tails). Never grows: green parcels (his approved/manual), parcels > 300,000 m2. Gaps > 20,000 m2 and
notches open to the outside (< 70 % of the gap outline against parcels) are not touched.
Receiving side: team CMP parcels first (already distorted), else the georeferenced sheet whose plots grow least in relative terms.
Writes fix_edgematch.json (+ edgematch_report.csv, edgematch_parcels.csv). Usage: python edge_match.py [x0 y0 x1 y1]"""
import sys, math, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd, geopandas as gpd
from shapely.geometry import Polygon, LineString, MultiLineString, Point, box
from shapely import make_valid
from shapely.ops import unary_union, split, linemerge
from shapely.strtree import STRtree

import os
LIVE = os.environ.get("CF_FILE", r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); PFX = os.environ.get("CF_PFX", "two")
P = gpd.read_file(LIVE).to_crs(32644)
if os.environ.get("EM_GEOM"):      # second pass: geometry from a preview file, attributes from the live layer (same row order)
    Q = gpd.read_file(os.environ["EM_GEOM"]).to_crs(32644); assert len(Q) == len(P) and (Q.kide.astype(str).values == P.kide.astype(str).values).all()
    P0 = P.copy(); P = P.set_geometry(Q.geometry.values)
P["kide"] = P.kide.astype(str); P["vc"] = P.village_code.astype(str); P["sv"] = P.kide.str.split("/").str[0]
WIN = box(*map(float, sys.argv[1:5])) if len(sys.argv) >= 5 else None
geom = {i: make_valid(g) for i, g in zip(P.index, P.geometry)}
orig = dict(geom)
live = {i: make_valid(g) for i, g in zip(P.index, (P0 if os.environ.get("EM_GEOM") else P).geometry)}
fixed = {i for i in P.index if str(P.status[i]).startswith("green") or geom[i].area > 300000}
isgeo = {i: str(P.source[i]).startswith("georef") for i in P.index}
U = make_valid(unary_union(list(geom.values())))
holes = [Polygon(r) for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
cl = U.buffer(8, join_style=2).buffer(-8, join_style=2)
strips = [p for p in getattr(make_valid(cl.difference(U)), "geoms", []) if p.geom_type == "Polygon"]
G = make_valid(unary_union([h for h in holes if h.area < 20000 and h.area / max(h.length / 2, 1e-9) <= 5.0] + strips))   # open areas wider than 5 m are outside the buffer: left
G = [g for g in getattr(G, "geoms", [G]) if g.geom_type == "Polygon" and g.area > 0.5 and g.area / max(g.length / 2, 1e-9) <= float(os.environ.get("EM_MAXW", 3.0))]   # seams only; layout roads / open land stay
if WIN is not None:
    G = [g for g in G if g.intersects(WIN)]
if os.environ.get("EM_GAPS"):     # externally prepared gap pieces (e.g. seam halves) instead of the layer's own holes
    G = [make_valid(g) for g in gpd.read_file(os.environ["EM_GAPS"]).geometry]
keys = list(geom)
T = STRtree([geom[k] for k in keys])
rep, skipped = [], []


def contact(i, g):
    return g.boundary.intersection(geom[i].buffer(0.05)).length


def lines_of(sh):
    if sh.geom_type == "LineString":
        return [sh]
    m = linemerge(sh) if sh.geom_type == "MultiLineString" else sh
    return [s for s in getattr(m, "geoms", [m]) if s.geom_type == "LineString"]


for g in sorted(G, key=lambda g: -g.area):
    nb = [keys[j] for j in T.query(g.buffer(0.1))]
    cont = {i: contact(i, g) for i in nb}
    cont = {i: c for i, c in cont.items() if c > 0.3}
    enclosed = sum(cont.values()) / g.length
    if enclosed < (0.5 if g.area < 100 else 0.7):
        skipped.append((round(g.area), round(g.centroid.x), round(g.centroid.y), f"open to the outside ({enclosed:.0%} against parcels)"))
        continue
    cand = [i for i in cont if i not in fixed]
    if not cand:
        skipped.append((round(g.area), round(g.centroid.x), round(g.centroid.y), "only approved/fixed parcels around"))
        continue
    groups = {}
    for i in cand:
        groups.setdefault((P.vc[i], P.sv[i]), []).append(i)

    width = g.area / max(g.length / 2, 1)

    def plan(recv):
        rays = []
        allp = list(cont)
        for a in allp:
            for b in allp:
                if a >= b:
                    continue
                sh = geom[a].boundary.intersection(geom[b].boundary)
                if sh.is_empty or sh.length < 0.05:
                    continue
                for s_ in lines_of(sh):
                    c = list(s_.coords)
                    for end, prev in ((c[0], c[1]), (c[-1], c[-2])):
                        if Point(end).distance(g) > 0.15:
                            continue
                        dx, dy = end[0] - prev[0], end[1] - prev[1]
                        L = math.hypot(dx, dy) or 1
                        far = 4 * width + 25
                        rays.append(LineString([(end[0] - dx / L * 0.2, end[1] - dy / L * 0.2), (end[0] + dx / L * far, end[1] + dy / L * far)]))
        pieces = [g]
        if rays:
            try:
                pieces = [p for p in split(g, MultiLineString(rays)).geoms if p.area > 0.01]
            except Exception:
                pieces = [g]
        give = {}
        for p in pieces:
            c = {i: contact(i, p) for i in recv}
            best = max(c, key=c.get)
            if c[best] < 0.05:
                oth = {i: contact(i, p) for i in cand}
                best = max(oth, key=oth.get)
            give.setdefault(best, []).append(p)
        growth = max(sum(q.area for q in v) / geom[i].area for i, v in give.items())
        # a piece much longer than the plot's own contact with the gap is a tail
        tail = max(max(q.length for q in v) / max(contact(i, g), 0.5) for i, v in give.items())
        return give, growth, tail, len(pieces)
    opts = []
    if os.environ.get("EM_RECV"):     # only these surveys may receive the gap
        groups = {k: v for k, v in groups.items() if str(k[1]) in os.environ["EM_RECV"].split(",")}
    for k, ids in sorted(groups.items(), key=lambda kv: -sum(cont[i] for i in kv[1])):
        give, growth, tail, npc = plan(ids)
        opts.append((growth <= float(os.environ.get("EM_MAXGROW", 0.30)) and tail <= 3.0, -sum(cont[i] for i in ids), growth, k, ids, give, npc))
    ok = [o for o in opts if o[0]]
    if not ok:
        skipped.append((round(g.area), round(g.centroid.x), round(g.centroid.y), "every side would grow > 30 % or get a tail - left open")); continue
    chosen = ok[0]
    _, _, growth, gk, recv, give, npc = chosen
    for i, ps in give.items():
        geom[i] = make_valid(unary_union([geom[i]] + [p.buffer(0.001, join_style=2) for p in ps]))
    pieces = [p for v in give.values() for p in v]
    if not ok:
        print(f"  gap {g.area:.0f} m2 at {g.centroid.x:.0f},{g.centroid.y:.0f}: every side grows > 30 % or makes a tail; used {gk} (max growth {growth:.0%})")
    rep.append({"gap_m2": round(g.area, 1), "x": round(g.centroid.x), "y": round(g.centroid.y), "width_m": round(width, 1),
                "receiving": f"{gk[0]} {gk[1]}" + (" (sheet)" if any(isgeo[i] for i in recv) else " (team CMP)"), "plots": len(recv), "pieces": len(pieces)})

changed = [i for i in keys if not geom[i].equals(orig[i])]
changed_live = lambda: [i for i in keys if not geom[i].equals(live[i])]
for i in changed:
    for j in T.query(geom[i]):
        j = keys[j]
        other = geom[j] if j in changed else orig[j]
        if j != i and geom[i].intersection(other).area > 1e-4 and (j not in changed or j < i):
            geom[i] = make_valid(geom[i].difference(other))
    g = geom[i]
    if g.geom_type != "Polygon":
        parts = [q for q in getattr(g, "geoms", [g]) if q.geom_type == "Polygon"]
        big = max(parts, key=lambda q: q.area)
        if g.area - big.area > 0.5:
            print(f"  {P.kide[i]}: {g.area - big.area:.1f} m2 in a detached part dropped")
        geom[i] = big
    geom[i] = Polygon(geom[i].exterior, [r for r in geom[i].interiors if Polygon(r).area > 0.5])

json.dump({"change": [], "add": [], "delete": []}, open(f"fix_{PFX}_edgematch.json", "w"))
R = pd.DataFrame(rep); R.to_csv(f"{PFX}_edgematch_report.csv", index=False)
ch = pd.DataFrame([{"kide": P.kide[i], "vc": P.vc[i], "status": str(P.status[i])[:14], "before": round(live[i].area, 1), "after": round(geom[i].area, 1),
                    "grow_m2": round(geom[i].area - live[i].area, 1), "grow_pct": round((geom[i].area / live[i].area - 1) * 100, 1)} for i in changed_live()])
ch = ch.sort_values("grow_pct", ascending=False) if len(ch) else ch
ch.to_csv(f"{PFX}_edgematch_parcels.csv", index=False)
print(f"gaps closed {len(R)} ({R.gap_m2.sum() if len(R) else 0:.0f} m2) | skipped {len(skipped)} | parcels changed {len(changed)}")
print("skipped:", skipped[:20])
print("largest relative growth:\n", ch.head(30).to_string(index=False))
out = {"change": [], "add": [], "delete": []}
for i in changed_live():
    g = geom[i]
    rp = live[i].representative_point()
    note = (str(P.merge_note[i] or "") + "; 2026-10-23: edge-matched onto its neighbour to close a gap").strip("; ")
    out["change"].append({"vc": P.vc[i], "old_kide": P.kide[i], "old_pt": list(rp.coords[0]),
                          "wkt": g.wkt,
                          "attrs": {"area_sqm": round(g.area, 1), "area_acre": round(g.area / 4046.856, 4), "merge_note": note}})
json.dump(out, open(f"fix_{PFX}_edgematch.json", "w"))
gpd.GeoDataFrame({"kide": [P.kide[i] for i in keys], "vc": [P.vc[i] for i in keys], "source": [P.source[i] for i in keys]},
                 geometry=[geom[i] for i in keys], crs=32644).to_file(f"{PFX}_edgematch_preview.gpkg", driver="GPKG")
U2 = make_valid(unary_union(list(geom.values())))
h2 = [Polygon(r).area for pp in getattr(U2, "geoms", [U2]) if pp.geom_type == "Polygon" for r in pp.interiors if Polygon(r).area < 20000]
cl2 = U2.buffer(8, join_style=2).buffer(-8, join_style=2)
s2 = [p for p in getattr(make_valid(cl2.difference(U2)), "geoms", []) if p.geom_type == "Polygon" and p.area > 0.5]
T2 = STRtree([geom[k] for k in keys])
ov = sum(geom[keys[a]].intersection(geom[keys[b]]).area for a in range(len(keys)) for b in T2.query(geom[keys[a]]) if b > a)
print(f"after: holes {len([h for h in h2 if h > 0.5])} ({sum(h2):.1f} m2), strips {len(s2)} ({sum(p.area for p in s2):.0f} m2), overlap total {ov:.2f} m2, invalid {sum(not geom[i].is_valid for i in changed)}")
