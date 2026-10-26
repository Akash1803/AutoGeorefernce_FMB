"""Close a seam strip by giving each side HALF of it: the strip is split at its midline (Voronoi equidistant line), each half is cut by
lines perpendicular to the strip axis at the points where the receiving plots meet, and every piece joins the plot it touches most.
Usage: half_assign.py <gaps.csv> <geojson> <row:surveysA/surveysB> <out_fix.json>"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, pandas as pd, numpy as np, shapely
from shapely import make_valid
from shapely.geometry import Polygon, LineString, Point, MultiPoint
from shapely.ops import unary_union, split, linemerge
from shapely.strtree import STRtree
CSV, F, SPEC, OUT = sys.argv[1:5]
P = gpd.read_file(F); P["kide"] = P.kide.astype(str); P["sv"] = P.kide.str.split("/").str[0]; P["vc"] = P.village_code.astype(str)
G = {i: make_valid(g) for i, g in zip(P.index, P.geometry)}; U = unary_union(list(G.values()))
holes = [Polygon(r) for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
D = pd.read_csv(CSV); r, sides = SPEC.split(":"); sa, sb = [set(s.split("+")) for s in sides.split("/")]; row = D.iloc[int(r)]
S = min(holes, key=lambda z: abs(z.area - row.area_m2) + z.distance(Point(row.x, row.y)))
near = P[P.intersects(S.buffer(0.05))]
side_of = {i: ("A" if near.sv[i] in sa else "B" if near.sv[i] in sb else None) for i in near.index}
UA = unary_union([G[i] for i in near.index if side_of[i] == "A"]); UB = unary_union([G[i] for i in near.index if side_of[i] == "B"])
ring = S.exterior; n = max(50, int(ring.length / 0.2))
pts = [ring.interpolate(t, normalized=True) for t in np.linspace(0, 1, n, endpoint=False)]
lab = ["A" if p.distance(UA) <= p.distance(UB) else "B" for p in pts]
cells = list(shapely.voronoi_polygons(MultiPoint(pts), extend_to=S.buffer(5)).geoms); T = STRtree(cells); own = {}
for k, p in enumerate(pts):
    c = T.query(p, predicate="within")
    if len(c): own.setdefault(int(c[0]), lab[k])
SA = make_valid(unary_union([cells[j] for j, l in own.items() if l == "A"]).intersection(S)); SB = make_valid(S.difference(SA))
mid = SA.boundary.intersection(SB.boundary); mid = linemerge(mid) if mid.geom_type != "LineString" else mid
mids = [m.simplify(0.03) for m in getattr(mid, "geoms", [mid]) if m.geom_type == "LineString" and m.length > 0.3]
def extend(m, d=3):
    c = np.array(m.coords); d0 = c[0] - c[1]; d1 = c[-1] - c[-2]
    return LineString([tuple(c[0] + d * d0 / np.linalg.norm(d0))] + [tuple(x) for x in c] + [tuple(c[-1] + d * d1 / np.linalg.norm(d1))])
pieces = [S]
for m in mids:
    pieces = [q for p in pieces for q in split(p, extend(m)).geoms if q.area > 0.005]
def contact(g, i): return g.buffer(0.02).intersection(G[i].boundary).length
halves = []
for p in pieces:
    ca = sum(contact(p, i) for i in near.index if side_of[i] == "A"); cb = sum(contact(p, i) for i in near.index if side_of[i] == "B")
    halves.append((p, "A" if ca >= cb else "B"))
print("strip %.1f m2 -> halves" % S.area, [(round(p.area, 1), s) for p, s in halves])
# strip axis for the perpendicular cuts
mr = shapely.minimum_rotated_rectangle(S); mc = np.array(mr.exterior.coords)[:4]; e = [mc[1] - mc[0], mc[2] - mc[1]]; axis = e[0] if np.linalg.norm(e[0]) > np.linalg.norm(e[1]) else e[1]
axis = axis / np.linalg.norm(axis); perp = np.array([-axis[1], axis[0]])
give = {}
for half, side in halves:
    recv = [i for i in near.index if side_of[i] == side and contact(half, i) > 0.05]
    # cut points: where two receiving plots meet on the half's boundary
    cuts = []
    for a in range(len(recv)):
        for b in range(a + 1, len(recv)):
            x = G[recv[a]].boundary.intersection(G[recv[b]].boundary)
            for q in getattr(x, "geoms", [x]):
                if q.is_empty: continue
                if q.geom_type == "Point" and q.distance(half) < 0.05: cuts.append(q)
                elif q.geom_type == "LineString":
                    for cc in (q.coords[0], q.coords[-1]):
                        if Point(cc).distance(half) < 0.05: cuts.append(Point(cc))
    parts = [half]
    for q in cuts:
        c = np.array(q.coords[0]); ln = LineString([tuple(c - 8 * perp), tuple(c + 8 * perp)])
        parts = [z for p in parts for z in split(p, ln).geoms if z.area > 0.005]
    for p in parts:
        cs = {i: contact(p, i) for i in recv}
        best = max(cs, key=cs.get) if cs and max(cs.values()) > 0.05 else None
        if best is None:                                   # no contact: join the neighbouring piece's plot
            nb = [(p.buffer(0.02).intersection(z).area, i2) for i2, zs in give.items() for z in zs]
            best = max(nb)[1] if nb else recv[0]
        give.setdefault(best, []).append(p)
out = {"change": [], "add": [], "delete": []}; rows = []
for i, ps in give.items():
    g = make_valid(unary_union([G[i]] + [p.buffer(0.001, join_style=2) for p in ps]))
    g = max([z for z in getattr(g, "geoms", [g]) if z.geom_type == "Polygon"], key=lambda z: z.area) if g.geom_type != "Polygon" else g
    g = Polygon(g.exterior, [h for h in g.interiors if Polygon(h).area > 0.5]).simplify(0.005)
    rows.append((P.kide[i], round(G[i].area, 1), round(g.area, 1), round((g.area / G[i].area - 1) * 100, 1)))
    note = (str(P.merge_note[i] if str(P.merge_note[i]) not in ("nan", "None") else "") + "; 2026-10-26: took its half of the seam strip to the neighbouring survey (edge matched to the midline)").strip("; ")
    out["change"].append({"vc": P.vc[i], "old_kide": P.kide[i], "old_pt": list(G[i].representative_point().coords[0]), "wkt": g.wkt,
                          "attrs": {"area_sqm": round(g.area, 1), "area_acre": round(g.area / 4046.856, 4), "merge_note": note}})
print("growth:", sorted(rows, key=lambda t: -t[3]))
json.dump(out, open(OUT, "w")); print("written", OUT, len(out["change"]), "| strip covered: %.1f of %.1f m2" % (unary_union([p for ps in give.values() for p in ps]).area, S.area))
