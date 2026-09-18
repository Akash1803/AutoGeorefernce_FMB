"""Village auto-georeferencing (spec 2026-09-17, revised): seeds = rigid poses recovered from the team's hand
files (pose only, never their geometry); mosaic chaining places the remaining sheets in the buffer; a whole-block
rigid adjustment settles all of them; a second chain pass picks up sheets that only match after the adjustment.
Puvi village vectors are a REFERENCE ONLY: the distance to the Puvi polygon of the same number is reported, it
never accepts, rejects or moves a placement. Writes <survey>_parcels_AUTOv2.gpkg (+ edges) and georef_status.csv.
Usage: python georef_village.py 35_04_077 [--seed-dir DIR] [--quiet]"""
import sys, os, glob, math, json, re, itertools, datetime, csv, argparse
import numpy as np, geopandas as gpd
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union
from shapely import affinity
from scipy.optimize import least_squares
from pyproj import Proj
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import autogeoref_demo as A
ap = argparse.ArgumentParser(); ap.add_argument("village"); ap.add_argument("--seed-dir", default=A.BK); ap.add_argument("--quiet", action="store_true")
args = ap.parse_args(); VILLAGE = args.village; A.set_village(VILLAGE); D = A.D; A.DEBUG = not args.quiet
ACRE = 0.000247105381; SIG_PAIR = 0.30; ADMIT = 1.0
key = lambda s: (int(re.match(r"\d+", s).group()), s)
def theta_of(R): return math.degrees(math.atan2(R[1, 0], R[0, 0]))
def R_of(theta): th = math.radians(theta); return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
# ---------------------------------------------------------------- 1. targets: buffer plots of this village with a sheet
buf = gpd.read_file("D:/Projects/Tambaram_Chengalpattu/Railway_Buffer_Vector_Plots.gpkg", layer="vector_in_buffer_30m")
in_buffer = sorted({str(s).strip() for s in buf[buf["village_code"].astype(str).str.strip() == VILLAGE]["survey_no"]}, key=key)
have_sheet = [sv for sv in in_buffer if os.path.exists(os.path.join(D, sv + "_parcels.geojson"))]
# railway land parcels (reference knowledge from the buffer plots: the track runs inside them) -> rail rule in A.rail_ok
_bp = buf[buf["village_code"].astype(str).str.strip() == VILLAGE].to_crs(32644)
A.RAIL_PARCELS = {str(r.survey_no).strip() for r in _bp.itertuples() if r.geometry.intersection(A.rail_line()).length >= A.RAIL_INSIDE_MIN}
print("railway land parcels (track runs inside, reference):", sorted(A.RAIL_PARCELS))
print("village %s: %d survey units in the 30 m buffer, %d with a sheet; missing sheets: %s" % (VILLAGE, len(in_buffer), len(have_sheet), [s for s in in_buffer if s not in have_sheet]))
# ---------------------------------------------------------------- 2. seeds: rigid pose from every hand file
seeds = {}
for sv in have_sheet:
    hb = sorted(glob.glob(os.path.join(args.seed_dir, sv + "_parcels_modified*.gpkg")), key=os.path.getmtime)
    if not hb: continue
    sheet, R, t, rms = A.seed_from_hand(sv, hb[-1])
    seeds[sv] = {"sheet": sheet, "R": R, "t": np.asarray(t, float), "theta": theta_of(R), "rms": rms, "file": os.path.basename(hb[-1])}
    print("seed %-4s from %-28s heading %7.2f  hand-vs-rigid rms %.2f m" % (sv, seeds[sv]["file"], seeds[sv]["theta"], rms))
sheets = {sv: seeds[sv]["sheet"] for sv in seeds}
pose = {sv: (seeds[sv]["theta"], seeds[sv]["t"]) for sv in seeds}          # current pose of every placed sheet
sig = {sv: (2.5, 3.0) for sv in seeds}                                    # prior sigmas: seeds strong, chained weak
notes = {sv: "" for sv in have_sheet}; chained_info = {}
# ---------------------------------------------------------------- 3. chain pass (used twice: before and after adjustment)
def chain_pass(pending, label):
    placed = {sv: (sheets[sv], R_of(pose[sv][0]), pose[sv][1]) for sv in pose}
    new = []
    while pending:
        best = None
        for sv in pending:
            r = A.place(sv, placed, 0.0)
            if r is None: continue
            rank = ((0 if r["single_edge"] else 1), r["conf"])       # multi-vertex or multi-anchor placements always go before single edges
            if best is None or rank > best[2]: best = (sv, r, rank)
        if best is None: break
        sv, r, _ = best; pending.remove(sv); new.append(sv)
        sheets[sv] = r["sheet"]; pose[sv] = (r["theta"], np.asarray(r["t"], float)); sig[sv] = (10.0, 10.0); chained_info[sv] = r
        if not r["single_edge"] or r["tentative_anchor"]:
            placed[sv] = (r["sheet"], r["R"], r["t"])
            if r["single_edge"]: A.TENTATIVE.add(sv); notes[sv] = "single edge, named by neighbour (tentative anchor)"
        else: notes[sv] = "single edge only, not used as anchor"
        if not r.get("supported", True): notes[sv] = (notes[sv] + " | " if notes[sv] else "") + "no printed neighbour number corroborates this chain"
        print("%s chained %-4s conf=%3d rms=%.2f m heading %7.2f anchors=%s rejected=%s %s" % (label, sv, r["conf"], r["rms"], r["theta"], [(a[0], a[1], round(a[2], 1)) for a in r["anchors"]], r["rejected"], notes.get(sv, "")))
    return new, pending
# ---------------------------------------------------------------- 4. whole-block rigid adjustment
def block_adjust():
    svs = sorted(pose, key=key); idx = {sv: i for i, sv in enumerate(svs)}
    outl = {sv: A.outline(sheets[sv]) for sv in svs}
    O = np.mean([pose[sv][1] for sv in svs], axis=0); x0 = np.zeros(3 * len(svs))
    for sv in svs: x0[3 * idx[sv]] = math.radians(pose[sv][0]); x0[3 * idx[sv] + 1:3 * idx[sv] + 3] = pose[sv][1] - O
    def T(sv, x, pt):
        th, tx, ty = x[3 * idx[sv]:3 * idx[sv] + 3]; c, s = math.cos(th), math.sin(th)
        return np.array([c * pt[0] - s * pt[1] + tx, s * pt[0] + c * pt[1] + ty])
    def apply_x(sv, x, geom):
        th, tx, ty = x[3 * idx[sv]:3 * idx[sv] + 3]
        return affinity.translate(affinity.rotate(geom, math.degrees(th), origin=(0, 0)), tx, ty)
    placed_geom = {sv: unary_union([apply_x(sv, x0, g) for _, g in sheets[sv]]).buffer(0) for sv in svs}
    start_c = {sv: np.array(placed_geom[sv].centroid.coords[0]) for sv in svs}
    sheet_c = {sv: np.array(unary_union([g for _, g in sheets[sv]]).centroid.coords[0]) for sv in svs}
    cands = []
    for a, b in itertools.combinations(svs, 2):
        if placed_geom[a].distance(placed_geom[b]) > 8.0: continue
        if A.adjacent_by_sheet(a, b) is False: continue          # transcribed sheets say these two are not neighbours
        best = None
        for L, pairs in A.common_chain(outl[a], [tuple(T(b, x0, p)) for p in outl[b]]):
            th, tx, ty = x0[3 * idx[b]:3 * idx[b] + 3]; c, s = math.cos(th), math.sin(th); loc = []
            for pa, pb in pairs:
                q = np.array(pb) - np.array([tx, ty]); loc.append((pa, (c * q[0] + s * q[1], -s * q[0] + c * q[1])))
            res = np.array([np.linalg.norm(T(a, x0, pa) - T(b, x0, pb)) for pa, pb in loc])
            if res.max() > (12.0 if len(loc) >= 3 else 6.0) or L < 8.0: continue     # hand poses can be 5-10 m off; multi-vertex chains are trusted further
            if best is None or L > best[0]: best = (L, len(loc), float(res.mean()), loc)
        if best: cands.append((a, b) + best)
    print("candidate boundaries:", [(a, b, round(L, 1), n, round(m, 2)) for a, b, L, n, m, _ in cands])
    obs = []
    def residuals(x):
        r = []
        for a, b, pa, pb in obs: r += list((T(a, x, pa) - T(b, x, pb)) / SIG_PAIR)
        for sv in svs:
            r += list((T(sv, x, sheet_c[sv]) - start_c[sv]) / sig[sv][0])
            r.append(((x[3 * idx[sv]] - x0[3 * idx[sv]] + math.pi) % (2 * math.pi) - math.pi) / math.radians(sig[sv][1]))
        return np.array(r)
    def solve(): return least_squares(residuals, x0, loss="huber", f_scale=1.5, x_scale="jac", xtol=1e-10, ftol=1e-10, gtol=1e-10)
    def chain_rms(x, c): return float(np.sqrt(np.mean([np.sum((T(c[0], x, pa) - T(c[1], x, pb)) ** 2) for pa, pb in c[5]])))
    used = [c for c in cands if c[3] >= 3]
    for c in used: obs += [(c[0], c[1], pa, pb) for pa, pb in c[5]]
    sol = solve(); x = sol.x; demoted = []
    while used:
        worst = max(used, key=lambda c: chain_rms(x, c))
        if chain_rms(x, worst) <= ADMIT: break
        demoted.append((worst[0], worst[1], round(chain_rms(x, worst), 2))); used.remove(worst)
        obs = [o for o in obs if (o[0], o[1]) != (worst[0], worst[1])]; sol = solve(); x = sol.x
    pending_c = [c for c in cands if c not in used]; rejected = []
    while pending_c:
        scored = sorted(((chain_rms(x, c), c) for c in pending_c), key=lambda t: t[0]); admitted = False
        for r0, c in scored:
            keep = obs; obs = obs + [(c[0], c[1], pa, pb) for pa, pb in c[5]]; sol_t = solve()
            wb = max([chain_rms(x, u) for u in used], default=0.0); wa = max([chain_rms(sol_t.x, u) for u in used], default=0.0)
            if chain_rms(sol_t.x, c) <= ADMIT and wa <= max(ADMIT, wb + 0.1):
                used.append(c); pending_c.remove(c); sol = sol_t; x = sol.x; admitted = True; break
            obs = keep
        if not admitted: rejected = [(c[0], c[1], round(chain_rms(x, c), 2)) for c in pending_c]; break
    print("boundaries used:", [(c[0], c[1], round(c[2], 1), c[3]) for c in used])
    print("boundaries demoted (sheets disagree):", demoted, "| chains rejected:", rejected)
    def overlap_total(xx):
        G = {sv: unary_union([apply_x(sv, xx, g) for _, g in sheets[sv]]).buffer(0) for sv in svs}
        return sum(G[a].intersection(G[b]).area for a, b in itertools.combinations(svs, 2))
    print("pairwise overlap between surveys: before adjustment %.1f m2 -> after %.1f m2" % (overlap_total(x0), overlap_total(x)))
    per = {}
    for a, b, pa, pb in obs: per.setdefault((a, b), []).append(float(np.linalg.norm(T(a, x, pa) - T(b, x, pb))))
    result = {}
    for sv in svs:
        th = math.degrees(x[3 * idx[sv]]); shift = float(np.linalg.norm(T(sv, x, sheet_c[sv]) - start_c[sv]))
        result[sv] = {"theta": th, "t": x[3 * idx[sv] + 1:3 * idx[sv] + 3] + O, "shift": shift, "dtheta": ((th - pose[sv][0]) + 180) % 360 - 180}
    return result, used, demoted, rejected, per, sol
# ---------------------------------------------------------------- 5. run: chain -> adjust -> chain again -> adjust
new1, waiting = chain_pass([sv for sv in have_sheet if sv not in seeds], "pass1")
adj, used, demoted, rejected, per, sol = block_adjust()
for sv, r in adj.items(): pose[sv] = (r["theta"], r["t"])
if waiting:
    new2, waiting = chain_pass(waiting, "pass2")
    if new2:
        for sv in new2: notes[sv] = (notes[sv] + " | " if notes[sv] else "") + "matched only after the block adjustment"
        adj, used, demoted, rejected, per, sol = block_adjust()
        for sv, r in adj.items(): pose[sv] = (r["theta"], r["t"])
if waiting: print("WAITING (no shared boundary with any placed sheet):", waiting)
# ---------------------------------------------------------------- 6. status, confidence, files
svs = sorted(pose, key=key); f = lambda v, full, zero: max(0.0, min(1.0, (zero - v) / (zero - full)))
proj = Proj("EPSG:32644"); run = datetime.datetime.now().isoformat(timespec="seconds"); status = []
for sv in svs:
    R = R_of(pose[sv][0]); t = pose[sv][1]; G = unary_union([A.apply(g, R, t) for _, g in sheets[sv]])
    mine = [(k, v) for k, v in per.items() if sv in k]; pr = [d for _, v in mine for d in v]
    nbrs = sorted({(k[1] if k[0] == sv else k[0]) for k, _ in mine}, key=key)
    L_tot = sum(c[2] for c in used if sv in (c[0], c[1])); n_pairs = sum(c[3] for c in used if sv in (c[0], c[1]))
    rms = float(np.sqrt(np.mean(np.square(pr)))) if pr else None
    conf = (40 * f(rms, 0.30, 1.5) + 30 * min(1.0, L_tot / 20.0) + 30 * min(2, len(nbrs)) / 2) if pr else 40.0
    if pr and n_pairs <= 2: conf = min(conf, 79)
    method = "seed" if sv in seeds else "chained"
    if sv in seeds: conf = max(conf, 80)
    pc = A.puvi_centroid(sv); puvi_m = round(float(G.centroid.distance(pc)), 1) if pc is not None else ""
    note = notes.get(sv, "")
    dis = [k for k in demoted + rejected if sv in (k[0], k[1])]
    if dis: note += (" | " if note else "") + "boundary disagreement with " + ",".join((k[1] if k[0] == sv else k[0]) + " (%.1f m)" % k[2] for k in dis)
    if not pr: note += (" | " if note else "") + "no shared boundary fitted; kept at " + ("hand" if sv in seeds else "chain") + " pose"
    if puvi_m != "" and puvi_m > 150: note += (" | " if note else "") + "far from Puvi reference (%.0f m) - check" % puvi_m; conf = min(conf, 60)
    row = {"survey": sv, "status": "seed" if sv in seeds else "placed", "method": method, "confidence": int(round(conf)), "boundary_rms_m": round(rms, 2) if rms is not None else "",
           "matched_length_m": round(L_tot, 1), "neighbours": ",".join(nbrs), "heading_deg": round(pose[sv][0], 2), "shift_from_start_m": round(adj[sv]["shift"], 2) if sv in adj else "",
           "hand_vs_rigid_rms_m": round(seeds[sv]["rms"], 2) if sv in seeds else "", "puvi_reference_m": puvi_m, "notes": note, "file": sv + "_parcels_AUTOv2.gpkg"}
    status.append(row); frs, ers = [], []
    for props, g in sheets[sv]:
        pg = A.apply(g, R, t); assert abs(pg.area - g.area) < 1e-6 and abs(pg.length - g.length) < 1e-6
        p = dict(props); a = pg.area
        p.update({"area_sqm": round(a, 3), "area_are": round(a / 100, 4), "area_hect": round(a / 1e4, 6), "area_acre": round(a * ACRE, 5), "area_cent": round(a * ACRE * 100, 3),
                  "fmb_area_sqm": round(a, 3), "fmb_area_acre": round(a * ACRE, 5), "fmb_area_cent": round(a * ACRE * 100, 3), "fmb_perimeter_m": round(pg.length, 3), "perimeter_m": round(pg.length, 3),
                  "georef_method": method + "+block_adjusted", "georef_confidence": row["confidence"], "georef_residual_m": rms, "anchored_to": row["neighbours"],
                  "heading_deg": row["heading_deg"], "puvi_reference_m": puvi_m, "georef_notes": note, "georef_run": run})
        frs.append({**p, "geometry": pg})
        c = list(pg.exterior.coords); c = c[::-1] if Polygon(c).exterior.is_ccw else c; st = max(range(len(c) - 1), key=lambda i: c[i][1]); c = c[st:-1] + c[:st]
        for i in range(len(c)):
            (x0_, y0_), (x1_, y1_) = c[i], c[(i + 1) % len(c)]; L = math.dist((x0_, y0_), (x1_, y1_)); grid = math.degrees(math.atan2(x1_ - x0_, y1_ - y0_)) % 360
            lon, lat = proj((x0_ + x1_) / 2, (y0_ + y1_) / 2, inverse=True); conv = proj.get_factors(lon, lat).meridian_convergence
            ers.append({"poly_id": props.get("poly_id"), "plot_no": props.get("plot_no"), "edge_no": i + 1, "length_m": round(L, 3), "bearing_grid_deg": round(grid, 2),
                        "bearing_true_deg": round((grid + conv) % 360, 2), "grid_convergence_deg": round(conv, 3), "geometry": LineString([(x0_, y0_), (x1_, y1_)])})
    out = os.path.join(D, sv + "_parcels_AUTOv2.gpkg")
    if os.path.exists(out): os.remove(out)
    gpd.GeoDataFrame(frs, geometry="geometry", crs="EPSG:32644").to_file(out, layer="parcels", driver="GPKG")
    gpd.GeoDataFrame(ers, geometry="geometry", crs="EPSG:32644").to_file(out, layer="edges", driver="GPKG")
for sv in waiting:
    status.append({"survey": sv, "status": "waiting", "method": "", "confidence": 0, "boundary_rms_m": "", "matched_length_m": 0, "neighbours": "", "heading_deg": "", "shift_from_start_m": "",
                   "hand_vs_rigid_rms_m": "", "puvi_reference_m": "", "notes": "no shared boundary with any placed sheet; seed by hand", "file": ""})
for f_ in glob.glob(os.path.join(D, "*_parcels_AUTOv2.gpkg")):
    if os.path.basename(f_).split("_")[0] not in svs: os.remove(f_)
with open(os.path.join(D, "georef_status.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(status[0].keys())); w.writeheader(); w.writerows(status)
print("\n%-5s %-7s %4s %7s %7s %-20s %8s %6s %7s %s" % ("svy", "status", "conf", "bnd_rms", "match_m", "neighbours", "heading", "shift", "puvi_m", "notes"))
for r in status: print("%-5s %-7s %4s %7s %7s %-20s %8s %6s %7s %s" % (r["survey"], r["status"], r["confidence"], r["boundary_rms_m"], r["matched_length_m"], r["neighbours"], r["heading_deg"], r["shift_from_start_m"], r["puvi_reference_m"], r["notes"]))
print("\nboundary residuals after adjustment:")
for (a, b), v in sorted(per.items(), key=lambda kv: -max(kv[1])): print("  %s-%s: %d pairs, rms %.2f m, max %.2f m" % (a, b, len(v), float(np.sqrt(np.mean(np.square(v)))), max(v)))
print("solver:", sol.status, "| nfev", sol.nfev, "| written", len(svs), "files to", D, "| waiting:", waiting)
