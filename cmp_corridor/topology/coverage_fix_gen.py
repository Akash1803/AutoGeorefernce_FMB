"""Exact coverage for CMP_my_Part without visible change (repeat until clean):
 round = node every boundary together on a 2e-8 deg grid -> polygonize -> faces back to their parcel:
   inside one parcel -> it | inside 2+ (hairline overlap) -> the bigger | inside none (hairline gap < 0.5 m2) -> cut at every point
   where the touching parcel changes, each piece to the parcel along it | gaps >= 0.5 m2 stay open (reported)
 the first round also snaps each parcel onto its neighbours within ~1.5 cm. Rounds repeat until the faces map 1:1.
-> fix_coverage.json, coverage_preview.gpkg"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, geopandas as gpd, shapely
from shapely import make_valid, set_precision
from shapely.geometry import Polygon, LineString, MultiLineString, Point
from shapely.ops import unary_union, polygonize, split
from shapely.strtree import STRtree
import os
F = os.environ.get("CF_FILE", r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); PFX = os.environ.get("CF_PFX", "two")
assert ("CF_FILE" in os.environ) == (PFX != "two"), "set CF_FILE and CF_PFX together"
GRID, MAXGAP, K, SNAP = 0.002, float(os.environ.get("CF_MAXGAP", 0.5)), 1.0, 0.015; WMIN = float(os.environ.get("CF_WMIN", 0.05))   # metres (EPSG:32644)
P = gpd.read_file(F); P["kide"] = P.kide.astype(str)
EMPTY = P[P.geometry.area < 0.01]
print("zero-area rows (to delete):", list(zip(EMPTY.village_code.astype(str), EMPTY.kide)))
EMPTY_DEL = [{"vc": str(v), "old_kide": k, "old_pt": list(g.centroid.coords[0])} for v, k, g in zip(EMPTY.village_code, EMPTY.kide, EMPTY.geometry)]
P = P[P.geometry.area >= 0.01].reset_index(drop=True)
def flat(g):
    out = []
    for z in getattr(g, "geoms", [g]):
        if z.geom_type == "Polygon": out.append(z)
        elif hasattr(z, "geoms"): out += flat(z)
    return out
G0, LOBES = [], []
for k_, g in zip(P.kide, P.geometry):
    if g.is_valid: G0.append(g); continue
    ps = sorted(flat(make_valid(g)), key=lambda z: -z.area); G0.append(ps[0])
    LOBES += [(k_, round(z.area, 2)) for z in ps[1:] if z.area > 0.01]
print("invalid parcels repaired: main body kept; lobes left out:", LOBES)
poly = lambda g: g if g.geom_type == "Polygon" else max([z for z in getattr(g, "geoms", [g]) if z.geom_type == "Polygon"], key=lambda z: z.area)

# snap onto neighbours (closes most hairlines before noding)
FROZEN = {i for i in range(len(P)) if str(P.village_code[i]) == "35_11_117" and P.kide[i] == "124"}
T0 = STRtree(G0); G = []
for i, g in enumerate(G0):
    if i in FROZEN: G.append(g); continue
    nb = [G0[int(j)] for j in T0.query(g.buffer(SNAP * 2)) if int(j) != i]
    s_ = poly(make_valid(shapely.snap(g, unary_union([q.boundary for q in nb]), SNAP))) if nb else g
    G.append(s_ if abs(s_.area - g.area) * K < 1.0 else g)


def cuts_for(f, nbi, Gc):
    cuts, ring = [], list(f.exterior.coords)[:-1]
    who = [tuple(sorted(j for j in nbi if Gc[j].boundary.distance(Point(p)) < GRID * 3)) for p in ring]
    for k, p in enumerate(ring):
        if len(who[k]) >= 2 or who[k] != who[k - 1]:
            a_, c_ = np.array(ring[k - 1]), np.array(ring[(k + 1) % len(ring)]); d = c_ - a_; n = np.hypot(*d)
            if n == 0: continue
            nx, ny = -d[1] / n, d[0] / n; L = 0.2
            cuts.append(LineString([(p[0] - nx * L, p[1] - ny * L), (p[0] + nx * L, p[1] + ny * L)]))
    return cuts


REAL_OV = []
def one_round(Gc):
    REAL_OV.clear()
    lines = unary_union(set_precision(unary_union([g.boundary for g in Gc]), GRID))
    faces = [f for f in polygonize(lines) if f.area > 0]
    T = STRtree(Gc); own = {i: [] for i in range(len(Gc))}; left = []; n_ov = n_gap = 0
    for f in faces:
        hit = [int(j) for j in T.query(f.representative_point(), predicate="within")]
        fz = [j for j in hit if j in FROZEN]
        if fz: own[fz[0]].append(f); continue          # a frozen parcel keeps all its ground
        if len(hit) == 1: own[hit[0]].append(f); continue
        if len(hit) > 1:
            wm = f.area / max(f.length / 2, 1e-12)
            if f.area < MAXGAP or wm < WMIN: own[max(hit, key=lambda j: Gc[j].area)].append(f); n_ov += 1
            else:
                for j in hit: own[j].append(f)
                REAL_OV.append((f.area, [P.kide[j] for j in hit], f.representative_point().coords[0]))
            continue
        width_m = f.area / max(f.length / 2, 1e-12)
        if f.area * K >= MAXGAP and width_m >= WMIN: left.append(f); continue      # real gap: kept open
        n_gap += 1
        nbi = [int(j) for j in T.query(f.buffer(GRID * 5))]
        pieces = [f]; cu = cuts_for(f, nbi, Gc)
        if cu:
            try: pieces = [q for q in split(f, MultiLineString(cu)).geoms if q.area > 0]
            except Exception: pass
        for q in pieces:
            c = [(q.boundary.intersection(Gc[j].boundary.buffer(GRID * 3)).length, j) for j in nbi if j not in FROZEN]
            c = [t for t in c if t[0] > 0]
            if c: own[max(c)[1]].append(q)
    new, bits = {}, []
    for i in range(len(Gc)):
        if not own[i]: new[i] = Gc[i]; print("  NO FACES for", P.kide[i]); continue
        u = unary_union(own[i])
        if u.geom_type == "MultiPolygon":
            parts = sorted(u.geoms, key=lambda q: -q.area); u = parts[0]; bits += [(i, q) for q in parts[1:]]
        new[i] = poly(make_valid(u))
    for i, q in bits:
        c = [(q.boundary.intersection(new[j].boundary.buffer(GRID * 3)).length, j) for j in new if j != i and new[j].distance(q) < GRID * 3]
        if c:
            j = max(c)[1]; new[j] = poly(make_valid(unary_union([new[j], q])))
    return [new[i] for i in range(len(Gc))], n_ov, n_gap, len(bits), left


for r in range(6):
    G, n_ov, n_gap, n_bits, left = one_round(G)
    print(f"round {r + 1}: hairline overlaps {n_ov}, hairline gaps {n_gap}, crumbs {n_bits}")
    if n_ov == 0 and n_gap == 0 and n_bits == 0: break


def report(Gc):
    T = STRtree(Gc)
    ov = sum(1 for a in range(len(Gc)) for b in T.query(Gc[a], predicate="overlaps") if b > a)
    U = unary_union(Gc); holes = sorted(Polygon(rr).area * K for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for rr in pp.interiors)
    close = sum(int((np.hypot(*np.diff(np.array(rr.coords), axis=0).T) < 0.001).sum()) for g in Gc for rr in [g.exterior] + list(g.interiors))
    inv = [P.kide[i] for i, g in enumerate(Gc) if not g.is_valid]
    return ov, holes, close, inv


ov, holes, close, inv = report(G)
Am = P; Bm = gpd.GeoSeries(G, crs=32644)
hd = [(round(a.hausdorff_distance(b), 3), str(P.kide[i]), round(b.area - a.area, 3)) for i, (a, b) in enumerate(zip(Am.geometry, Bm))]
hd.sort(reverse=True)
print(f"after: overlapping pairs {ov} | holes {len(holes)} (largest {[round(h, 2) for h in holes[-5:]]}, tiny < 0.5 m2: {sum(h < 0.5 for h in holes)}) | close vertices {close} | invalid {inv[:5]}")
print("largest edge shifts (m, kide, area change m2):", hd[:6])
print("open gaps kept:", [round(f.area * K, 1) for f in sorted(left, key=lambda f: -f.area)][:6])
out = {"change": [], "add": [], "delete": []}
for i, g in enumerate(G):
    if g.equals_exact(G0[i], 0): continue
    a = Bm.iloc[i].area
    out["change"].append({"vc": str(P.village_code[i]), "old_kide": P.kide[i], "old_pt": list(G0[i].representative_point().coords[0]), "wkt": g.wkt,
                          "attrs": {"area_sqm": round(a, 1), "area_acre": round(a / 4046.856, 4)}})
out["delete"] = EMPTY_DEL
json.dump(out, open(f"fix_{PFX}_coverage.json", "w")); print("prepared", len(out["change"]))
gpd.GeoDataFrame({"kide": P.kide, "village_code": P.village_code.astype(str)}, geometry=G, crs=32644).to_file(f"{PFX}_coverage_preview.gpkg", driver="GPKG")
print("REAL overlaps kept (m2, parcels):", sorted([(round(a, 1), k) for a, k, c in REAL_OV], reverse=True)[:80])
json.dump([(a, k, c) for a, k, c in REAL_OV], open(f"{PFX}_real_overlaps.json", "w"))
