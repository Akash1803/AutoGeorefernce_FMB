"""Whole-block rigid least-squares adjustment (approach C of the design) for one village block.
Inputs : sheet polygons (local metres), manual poses recovered by manual_poses.py (weak prior),
         shared-boundary vertex pairs found between neighbouring sheets (strong observation).
Output : <survey>_parcels_AUTOv2.gpkg (rigid, scale 1) + block_adjust_report.csv
Usage  : python block_adjust.py 35_04_052"""
import sys, os, json, math, pickle, itertools, datetime, csv, numpy as np, geopandas as gpd
from shapely.geometry import shape, Polygon, LineString, Point
from shapely.ops import unary_union
from shapely import affinity
from scipy.optimize import least_squares
from pyproj import Proj
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import autogeoref_demo as A
VILLAGE = sys.argv[1]; A.set_village(VILLAGE); D = A.D
poses = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "manual_poses_%s.pkl" % VILLAGE), "rb"))
SIG_PAIR, SIG_POS, SIG_HEAD = 0.30, 2.5, 3.0          # m, m, degrees
ACRE = 0.000247105381
svs = sorted(poses); idx = {sv: i for i, sv in enumerate(svs)}
sheets = {sv: A.load_sheet(sv) for sv in svs}
outl = {sv: A.outline(sheets[sv]) for sv in svs}
def T(sv, x, pt):
    th, tx, ty = x[3 * idx[sv]:3 * idx[sv] + 3]; c, s = math.cos(th), math.sin(th)
    return np.array([c * pt[0] - s * pt[1] + tx, s * pt[0] + c * pt[1] + ty])
def apply_x(sv, x, geom):
    th, tx, ty = x[3 * idx[sv]:3 * idx[sv] + 3]
    return affinity.translate(affinity.rotate(geom, math.degrees(th), origin=(0, 0)), tx, ty)
x0 = np.zeros(3 * len(svs))
for sv in svs:
    x0[3 * idx[sv]] = math.radians(poses[sv]["theta"]); x0[3 * idx[sv] + 1:3 * idx[sv] + 3] = poses[sv]["t"]
O = np.mean([poses[sv]["t"] for sv in svs], axis=0)          # local origin: keeps the unknowns small for the solver
for sv in svs: x0[3 * idx[sv] + 1:3 * idx[sv] + 3] -= O
# 1. shared-boundary pairs: for every two surveys whose manual placements come within 5 m, find chains whose
#    residual under the manual poses is small (they are the real shared boundaries)
placed_geom = {sv: unary_union([apply_x(sv, x0, g) for _, g in sheets[sv]]) for sv in svs}
manual_centroid = {sv: np.array(placed_geom[sv].centroid.coords[0]) for sv in svs}
sheet_centroid = {sv: np.array(unary_union([g for _, g in sheets[sv]]).centroid.coords[0]) for sv in svs}
cands = []   # (a, b, L, n_pairs, mean_res_manual, loc)
for a, b in itertools.combinations(svs, 2):
    dist = placed_geom[a].distance(placed_geom[b])
    if dist > 5.0: continue
    best = None; best_any = None
    for L, pairs in A.common_chain(outl[a], [tuple(T(b, x0, p)) for p in outl[b]]):
        # pairs: (a local vertex, b PLACED vertex) -> recover b local vertex by inverse transform
        th, tx, ty = x0[3 * idx[b]:3 * idx[b] + 3]; c, s = math.cos(th), math.sin(th)
        loc = []
        for pa, pb in pairs:
            q = np.array(pb) - np.array([tx, ty]); loc.append((pa, (c * q[0] + s * q[1], -s * q[0] + c * q[1])))
        res = np.array([np.linalg.norm(T(a, x0, pa) - T(b, x0, pb)) for pa, pb in loc])
        if best_any is None or L > best_any[0]: best_any = (round(L, 1), len(loc), round(float(res.max()), 2))
        if res.max() > 6.0 or L < 8.0: continue
        cand = (L, len(loc), float(res.mean()), loc)
        if best is None or (cand[0] > best[0]): best = cand
    if best: cands.append((a, b) + best)
    else: print("   no usable chain for %s-%s (manual distance %.1f m; longest equal-length chain (L, pairs, max res) = %s)" % (a, b, dist, best_any))
print("candidate boundaries:", [(a, b, round(L, 1), n, round(m, 2)) for a, b, L, n, m, _ in cands])
obs = []
def residuals(x):
    r = []
    for a, b, pa, pb in obs: r += list((T(a, x, pa) - T(b, x, pb)) / SIG_PAIR)
    for sv in svs:
        r += list((T(sv, x, sheet_centroid[sv]) - manual_centroid[sv]) / SIG_POS)
        r.append(((x[3 * idx[sv]] - x0[3 * idx[sv]] + math.pi) % (2 * math.pi) - math.pi) / math.radians(SIG_HEAD))
    return np.array(r)
def solve(): return least_squares(residuals, x0, loss="huber", f_scale=1.5, x_scale="jac", xtol=1e-10, ftol=1e-10, gtol=1e-10)
def chain_rms(x, c): return float(np.sqrt(np.mean([np.sum((T(c[0], x, pa) - T(c[1], x, pb)) ** 2) for pa, pb in c[5]])))
# trusted first: chains of >= 3 vertex pairs (two or more consecutive equal edges cannot slide or mirror)
ADMIT_DEMOTE = 1.0
used = [c for c in cands if c[3] >= 3]
for c in used: obs += [(c[0], c[1], pa, pb) for pa, pb in c[5]]
sol = solve(); x = sol.x
# demote trusted chains that the block cannot satisfy (the two sheets disagree on the corner angle): worst first
demoted = []
while used:
    worst = max(used, key=lambda c: chain_rms(x, c))
    if chain_rms(x, worst) <= ADMIT_DEMOTE: break
    demoted.append((worst[0], worst[1], round(chain_rms(x, worst), 2))); used.remove(worst)
    obs = [o for o in obs if (o[0], o[1]) != (worst[0], worst[1])]; sol = solve(); x = sol.x
print("multi-vertex chains demoted (block cannot fit them within %.1f m):" % ADMIT_DEMOTE, demoted)
# then admit single-edge chains one at a time, best-agreeing first, only while they fit the block within ADMIT_M
ADMIT_M = 1.0; pending = [c for c in cands if c[3] < 3] + [c for c in cands if (c[0], c[1], round(chain_rms(x, c), 2)) in demoted]; rejected = []
while pending:
    scored = sorted(((chain_rms(x, c), c) for c in pending), key=lambda t: t[0]); admitted = False
    for r0, c in scored:                                   # trial: add, re-solve, keep only if it and the others still fit
        trial = obs + [(c[0], c[1], pa, pb) for pa, pb in c[5]]; keep = obs; obs = trial; sol_t = solve()
        worst_before = max([chain_rms(x, u) for u in used], default=0.0)
        worst_used = max([chain_rms(sol_t.x, u) for u in used], default=0.0)
        if chain_rms(sol_t.x, c) <= ADMIT_M and worst_used <= max(ADMIT_M, worst_before + 0.1):
            used.append(c); pending.remove(c); sol = sol_t; x = sol.x; admitted = True; break
        obs = keep
    if not admitted: rejected = [(c[0], c[1], round(chain_rms(x, c), 2)) for c in pending]; break
print("boundaries used:", [(c[0], c[1], round(c[2], 1), c[3]) for c in used])
print("single-edge chains rejected (disagree with the block by > %.1f m rms; sheet inconsistency, review):" % ADMIT_M, rejected)
# manual vs adjusted: total pairwise overlap area between surveys
def overlap_total(xx):
    G = {sv: unary_union([apply_x(sv, xx, g) for _, g in sheets[sv]]).buffer(0) for sv in svs}
    return sum(G[a].intersection(G[b]).area for a, b in itertools.combinations(svs, 2))
print("pairwise overlap area between surveys: manual %.1f m2 -> adjusted %.1f m2" % (overlap_total(x0), overlap_total(x)))
# 2. report
rows = []
for sv in svs:
    th = math.degrees(x[3 * idx[sv]]); dth = ((th - poses[sv]["theta"]) + 180) % 360 - 180
    shift = float(np.linalg.norm(T(sv, x, sheet_centroid[sv]) - manual_centroid[sv]))
    pr = [np.linalg.norm(T(a, x, pa) - T(b, x, pb)) for a, b, pa, pb in obs if sv in (a, b)]
    rows.append({"survey": sv, "heading_manual": round(poses[sv]["theta"], 2), "heading_adjusted": round(th, 2), "heading_change_deg": round(dth, 2),
                 "shift_from_manual_m": round(shift, 2), "manual_rigid_rms_m": round(poses[sv]["rms"], 2), "manual_scale": round(poses[sv]["scale"], 4),
                 "boundary_pairs": len(pr), "boundary_rms_m": round(float(np.sqrt(np.mean(np.square(pr)))), 2) if pr else None, "boundary_max_m": round(float(max(pr)), 2) if pr else None,
                 "neighbours_used": ",".join(sorted({(b if a == sv else a) for a, b, _, _ in obs if sv in (a, b)}))})
# pair residuals after adjustment (sheet inconsistencies show up here)
after = {}
for a, b, pa, pb in obs: after.setdefault((a, b), []).append(np.linalg.norm(T(a, x, pa) - T(b, x, pb)))
pair_rows = [{"a": a, "b": b, "pairs": len(v), "rms_m": round(float(np.sqrt(np.mean(np.square(v)))), 2), "max_m": round(float(max(v)), 2)} for (a, b), v in after.items()]
# 3. write v2 files
proj = Proj("EPSG:32644"); run = datetime.datetime.now().isoformat(timespec="seconds")
for sv in svs:
    prow = next(r for r in rows if r["survey"] == sv); frs, ers = [], []
    for props, g in sheets[sv]:
        pg = affinity.translate(apply_x(sv, x, g), O[0], O[1]); assert abs(pg.area - g.area) < 1e-6 and abs(pg.length - g.length) < 1e-6
        p = dict(props); a = pg.area
        p.update({"area_sqm": round(a, 3), "area_are": round(a / 100, 4), "area_hect": round(a / 1e4, 6), "area_acre": round(a * ACRE, 5), "area_cent": round(a * ACRE * 100, 3),
                  "fmb_area_sqm": round(a, 3), "fmb_area_acre": round(a * ACRE, 5), "fmb_area_cent": round(a * ACRE * 100, 3), "fmb_perimeter_m": round(pg.length, 3), "perimeter_m": round(pg.length, 3),
                  "georef_method": "block_adjusted", "georef_confidence": int(max(0, min(100, 100 - 20 * (prow["boundary_rms_m"] or 0) / SIG_PAIR - 5 * prow["shift_from_manual_m"]))) if prow["boundary_pairs"] else 40,
                  "georef_residual_m": prow["boundary_rms_m"], "anchored_to": prow["neighbours_used"], "shift_from_manual_m": prow["shift_from_manual_m"], "heading_deg": prow["heading_adjusted"], "georef_run": run})
        frs.append({**p, "geometry": pg})
        c = list(pg.exterior.coords); c = c[::-1] if Polygon(c).exterior.is_ccw else c; st = max(range(len(c) - 1), key=lambda i: c[i][1]); c = c[st:-1] + c[:st]
        for i in range(len(c)):
            (x0_, y0_), (x1_, y1_) = c[i], c[(i + 1) % len(c)]; L = math.dist((x0_, y0_), (x1_, y1_)); grid = math.degrees(math.atan2(x1_ - x0_, y1_ - y0_)) % 360
            lon, lat = proj((x0_ + x1_) / 2, (y0_ + y1_) / 2, inverse=True); conv = proj.get_factors(lon, lat).meridian_convergence
            ers.append({"poly_id": props.get("poly_id"), "plot_no": props.get("plot_no"), "edge_no": i + 1, "length_m": round(L, 3), "bearing_grid_deg": round(grid, 2), "bearing_true_deg": round((grid + conv) % 360, 2), "geometry": LineString([(x0_, y0_), (x1_, y1_)])})
    out = os.path.join(D, sv + "_parcels_AUTOv2.gpkg")
    if os.path.exists(out): os.remove(out)
    gpd.GeoDataFrame(frs, geometry="geometry", crs="EPSG:32644").to_file(out, layer="parcels", driver="GPKG")
    gpd.GeoDataFrame(ers, geometry="geometry", crs="EPSG:32644").to_file(out, layer="edges", driver="GPKG")
with open(os.path.join(D, "block_adjust_report.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\n%-5s %8s %8s %7s %8s %9s %6s %8s %8s  %s" % ("svy", "head_man", "head_adj", "d_head", "shift_m", "man_rms", "pairs", "bnd_rms", "bnd_max", "neighbours"))
for r in rows: print("%-5s %8.2f %8.2f %7.2f %8.2f %9.2f %6d %8s %8s  %s" % (r["survey"], r["heading_manual"], r["heading_adjusted"], r["heading_change_deg"], r["shift_from_manual_m"], r["manual_rigid_rms_m"], r["boundary_pairs"], r["boundary_rms_m"], r["boundary_max_m"], r["neighbours_used"]))
print("\nshared-boundary residuals after adjustment (sheet-vs-sheet disagreement where large):")
for r in sorted(pair_rows, key=lambda r: -r["max_m"]): print("  %s-%s: %d pairs, rms %.2f m, max %.2f m" % (r["a"], r["b"], r["pairs"], r["rms_m"], r["max_m"]))
print("\nsolver:", sol.status, sol.message, "| cost", round(float(sol.cost), 1), "| nfev", sol.nfev, "| optimality", round(float(sol.optimality), 4))
