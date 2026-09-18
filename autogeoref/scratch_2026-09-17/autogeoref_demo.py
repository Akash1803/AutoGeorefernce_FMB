"""PROTOTYPE / DEMO of the mosaic-chaining georeferencing (design spec 2026-09-17).
Seeds with one hand-placed sheet (rigidly re-fitted so no stretch is inherited), then places
neighbours by matching shared boundary edge lengths and a rigid (rotate + shift) fit.
Writes <survey>_parcels_AUTODEMO.gpkg files only - never touches the team's _modified files."""
import os, sys, json, math, glob, datetime
import numpy as np, geopandas as gpd, pandas as pd
from shapely.geometry import shape, Polygon, LineString, Point
from shapely.ops import unary_union
from shapely import affinity
from pyproj import Proj

VILLAGE = "35_04_077"
D = "D:/Projects/Tambaram_Chengalpattu/FMB_Vector/" + VILLAGE
SEED_FILE = None
def set_village(code):
    global VILLAGE, D
    VILLAGE = code; D = "D:/Projects/Tambaram_Chengalpattu/FMB_Vector/" + code
BK = r"D:\Projects\Tambaram_Chengalpattu\_logs\backup_georef_20260917"
TOL_ABS, TOL_REL = 0.30, 0.01
ACRE = 0.000247105381

def load_sheet(sv):
    g = json.load(open(os.path.join(D, sv + "_parcels.geojson"), encoding="utf-8"))
    return [(f["properties"], shape(f["geometry"])) for f in g["features"]]

def outline(polys):
    """exterior ring of the dissolved sketch with collinear vertices merged (angle < 1 deg)."""
    U = unary_union([p for _, p in polys]).buffer(0)
    U = max(U.geoms, key=lambda q: q.area) if U.geom_type == "MultiPolygon" else U
    c = list(U.exterior.coords)[:-1]
    if Polygon(c).exterior.is_ccw: c = c[::-1]                    # clockwise
    out = []
    n = len(c)
    for i in range(n):
        a, b, d = c[i - 1], c[i], c[(i + 1) % n]
        ang = abs((math.degrees(math.atan2(d[1] - b[1], d[0] - b[0]) - math.atan2(b[1] - a[1], b[0] - a[0])) + 180) % 360 - 180)
        if ang > 0.1: out.append(b)          # lesson: no collinear merge (only numerical duplicates dropped)
    return out

def edges(v):
    return [math.dist(v[i], v[(i + 1) % len(v)]) for i in range(len(v))]

def turn_angles(v):
    """turn angle (deg) at each vertex of a closed outline, rotation-invariant."""
    n = len(v); out = []
    for i in range(n):
        a, b, d = v[i - 1], v[i], v[(i + 1) % n]
        out.append(abs((math.degrees(math.atan2(d[1] - b[1], d[0] - b[0]) - math.atan2(b[1] - a[1], b[0] - a[0])) + 180) % 360 - 180))
    return out

PASS_DEG = 5.0
BOTH_DIRS = True
MIN_CHAIN = 8.0
DEBUG = True      # a corner turning less than this is treated as straight (drawn slightly differently on the two sheets)

def common_chain(vq, vp):
    """Longest shared boundary between outline q and (placed) outline p.
    Walk both boundaries from a start corner pair by running distance. Corners at the same
    running distance on both sheets are matched; a corner on one sheet that falls inside a
    straight stretch of the other is passed through if its turn is < PASS_DEG, otherwise the
    boundaries diverge there and the walk ends. Returns (matched_length_m, pairs) or (0, None)."""
    def tol(x): return TOL_ABS + TOL_REL * max(x, 1.0)
    eq = edges(vq); nq = len(eq); tq = turn_angles(vq); best = (0.0, None); results = []; seen = set()
    for rev in ((True, False) if BOTH_DIRS else (True,)):
        vp2 = vp[::-1] if rev else vp; ep = edges(vp2); np_ = len(ep); tp = turn_angles(vp2)
        for i in range(nq):
            for j in range(np_):
                pairs = [(vq[i], vp2[j])]; cq = cp = 0.0; iq = ip = 0; matched = 0.0
                while iq < nq and ip < np_:
                    nq_next, np_next = cq + eq[(i + iq) % nq], cp + ep[(j + ip) % np_]
                    if abs(nq_next - np_next) <= tol(max(nq_next, np_next)):
                        cq, cp, iq, ip = nq_next, np_next, iq + 1, ip + 1
                        pairs.append((vq[(i + iq) % nq], vp2[(j + ip) % np_])); matched = min(cq, cp)
                        if iq >= nq or ip >= np_: break
                    elif nq_next < np_next:
                        iq += 1; cq = nq_next
                        if tq[(i + iq) % nq] > PASS_DEG: break        # q really turns here, p is straight -> diverge
                    else:
                        ip += 1; cp = np_next
                        if tp[(j + ip) % np_] > PASS_DEG: break
                if len(pairs) >= 2 and matched >= MIN_CHAIN:
                    key = tuple(sorted((round(a[0], 2), round(a[1], 2), round(b[0], 2), round(b[1], 2)) for a, b in pairs))
                    if key not in seen: seen.add(key); results.append((matched, pairs))
    results.sort(key=lambda r: -r[0])
    return results

def rigid_fit(pairs):
    """rotation + shift with scale fixed at 1 (Procrustes without scale)."""
    P = np.array([p for p, _ in pairs]); Q = np.array([q for _, q in pairs])
    pc, qc = P.mean(0), Q.mean(0); H = (P - pc).T @ (Q - qc); U, S, Vt = np.linalg.svd(H); R = Vt.T @ U.T
    if np.linalg.det(R) < 0: Vt[-1] *= -1; R = Vt.T @ U.T
    t = qc - R @ pc; res = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
    return R, t, float(np.sqrt((res ** 2).mean())), float(res.max())

def apply(geom, R, t):
    th = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    return affinity.translate(affinity.rotate(geom, th, origin=(0, 0)), t[0], t[1])

def seed_from_hand(sv, hand_file=None):
    """rigid re-fit of the sheet onto the team's hand placement (drops the hand stretch)."""
    sheet = load_sheet(sv)
    if hand_file is None:
        hand_file = glob.glob(os.path.join(BK, sv + "_parcels_modified*.gpkg"))[0]
    hand = gpd.read_file(hand_file)
    pairs = []
    if "poly_id" in hand.columns:
        byid = {int(r.poly_id): r.geometry for r in hand.itertuples() if r.poly_id is not None}
        for props, g in sheet:
            gb = byid.get(int(props.get("poly_id", -1)))
            if gb is not None:
                ca, cb = list(g.exterior.coords)[:-1], list(gb.exterior.coords)[:-1]
                if len(ca) == len(cb): pairs += list(zip(ca, cb))
    if len(pairs) < 3:
        A = sorted([p for _, p in sheet], key=lambda g: -g.area); B = sorted(hand.geometry.tolist(), key=lambda g: -g.area)
        for ga, gb in zip(A, B):
            ca, cb = list(ga.exterior.coords)[:-1], list(gb.exterior.coords)[:-1]
            if len(ca) == len(cb): pairs += list(zip(ca, cb))
    R, t, rms, mx = rigid_fit(pairs)
    return sheet, R, t, rms

NB = {}
PUVI = {}
TENTATIVE = set()
def puvi_centroid(sv):
    """centroid of the village vector polygon(s) with this survey number (exact, else numeric part); None if absent."""
    import re, glob as _g
    if VILLAGE not in PUVI:
        d_, t, v = VILLAGE.split("_"); shp = _g.glob(os.path.join("D:/Data/Puvi_data/vector", d_, t, v, "vector", "*.shp"))
        pv = gpd.read_file(shp[0]).to_crs(32644) if shp else None
        d = {}
        if pv is not None:
            pv["key"] = pv["survey_no"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
            for k, grp in pv.groupby("key"): d[k] = unary_union(grp.geometry.tolist()).centroid
        PUVI[VILLAGE] = d
    d = PUVI[VILLAGE]; k = sv.upper()
    if k in d: return d[k]
    num = re.match(r"^\d+", k)
    return d.get(num.group(0)) if num else None
HEADING_GATE_DEG = 360.0   # lesson (Thailavaram): sheets of one village have no common heading; gate disabled
PUVI_GATE_M = None      # user rule: Puvi is a reference only, never used to accept or reject a placement
NBPOS = {}
LABEL_REACH = 80.0   # a neighbour's number is written within this distance (m) of the shared boundary
SIDE_VEC = {"N": (0, 1), "NE": (1, 1), "E": (1, 0), "SE": (1, -1), "S": (0, -1), "SW": (-1, -1), "W": (-1, 0), "NW": (-1, 1)}
NB_OVERRIDE = {}
def _override():
    if VILLAGE not in NB_OVERRIDE:
        f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nb_override_%s.json" % VILLAGE)
        NB_OVERRIDE[VILLAGE] = json.load(open(f, encoding="utf-8")) if os.path.exists(f) else None
    return NB_OVERRIDE[VILLAGE]
def _norm(num):
    import re as _re
    n = str(num).upper().replace(" ", "").replace("S.NO", "").replace("S.NO.", "").split("/")[0]
    return n if _re.match(r"^\d+[A-Z]?$", n) else None
def sheet_neighbours(sv):
    """neighbour survey numbers printed on the sheet: glyph reader, merged with the independent transcription if present."""
    if sv not in NB:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import neighbour_labels as NL
        pdf = os.path.join("D:/Projects/Tambaram_Chengalpattu/FMB_Sketches/" + VILLAGE, sv + ".pdf")
        labs, U, polys, M, H = NL.neighbour_labels(pdf)
        glyph = {l["text"].upper(): (l["pt"][0] * M, l["pt"][1] * M) for l in labs if not l["inside"] and l["dist_m"] < 60}
        ov = _override()
        if ov is None or sv not in ov:
            NB[sv] = set(glyph); NBPOS[sv] = glyph
        else:
            read = {}
            for e in ov[sv]:
                n = _norm(e.get("number", ""))
                if n and (e.get("readers", 1) >= 2 or e.get("confidence", "high") == "high"): read[n] = e.get("side", "")
            # glyph fragments that are a prefix of a transcribed number (17 -> 170) are replaced; other glyph reads kept
            keep = {g: pos for g, pos in glyph.items() if not any(r != g and r.startswith(g) for r in read)}
            U_ = unary_union([p for _, p in load_sheet(sv)]); c = U_.centroid; minx, miny, maxx, maxy = U_.bounds
            pos = dict(keep)
            for n, side in read.items():
                if n in pos: continue
                dx, dy = SIDE_VEC.get(side.upper(), (0, 0)); L_ = math.hypot(dx, dy) or 1.0
                pos[n] = (c.x + dx / L_ * ((maxx - minx) / 2 + 10), c.y + dy / L_ * ((maxy - miny) / 2 + 10))
            NB[sv] = set(pos); NBPOS[sv] = pos
    return NB[sv]

RAIL = None; RAIL_PARCELS = set(); RAIL_CROSS_MAX = 5.0; RAIL_INSIDE_MIN = 50.0
def rail_line():
    global RAIL
    if RAIL is None:
        r = gpd.read_file("D:/Projects/Tambaram_Chengalpattu/Tambaram_Chengalpattu_Railway.gpkg", layer="rail_line").to_crs(32644)
        RAIL = unary_union(list(r.geometry))
    return RAIL
def rail_ok(sv, geom):
    """soft rule: railway land should contain the track; other parcels are penalised per metre of track inside them and
    vetoed only when the track runs through more than half of their longest dimension."""
    cross = geom.intersection(rail_line()).length
    if sv in RAIL_PARCELS: return True, cross
    c = list(geom.minimum_rotated_rectangle.exterior.coords); longest = max(math.dist(c[i], c[i + 1]) for i in range(len(c) - 1))
    return True, cross                      # score-only: strips beside the line often predate the railway land and lie under today's tracks

def label_conflicts(sv, sheet, R, t, placed):
    """(name, which, distance) for every printed number linking sv with a placed sheet that lands > LABEL_REACH away."""
    q_geom = unary_union([apply(p, R, t) for _, p in sheet]); out = []
    for name, (pp, Rp, tp) in placed.items():
        if name.upper() in NBPOS.get(sv, {}):
            d = Point(apply(Point(NBPOS[sv][name.upper()]), R, t)).distance(unary_union([apply(p, Rp, tp) for _, p in pp]))
            if d > LABEL_REACH: out.append((name, "q-label", round(d)))
        if sv.upper() in NBPOS.get(name, {}):
            d = Point(apply(Point(NBPOS[name][sv.upper()]), Rp, tp)).distance(q_geom)
            if d > LABEL_REACH: out.append((name, "p-label", round(d)))
    return out

def adjacent_by_sheet(a, b):
    ov = _override()
    if ov is None or a not in ov or b not in ov: return None
    return b.upper() in sheet_neighbours(a) or a.upper() in sheet_neighbours(b)

def place(sv, placed, prior_deg):
    sheet = load_sheet(sv); vq = outline(sheet); cands = []
    mine_nb = sheet_neighbours(sv)
    _rg = [unary_union([apply(p, Rp, tp) for _, p in pp]) for nm, (pp, Rp, tp) in placed.items() if nm in RAIL_PARCELS]
    rail_u = unary_union(_rg) if _rg else None
    def ov_area(geom_a, geom_b):
        inter = geom_a.intersection(geom_b)
        if rail_u is not None and not inter.is_empty: inter = inter.difference(rail_u)
        return inter.area
    for name, (polys_p, Rp, tp) in placed.items():
        vp = outline([(None, apply(p, Rp, tp)) for _, p in polys_p]); anchor_geom = unary_union([apply(p, Rp, tp) for _, p in polys_p])
        labelled = name.upper() in mine_nb or sv.upper() in sheet_neighbours(name)
        if adjacent_by_sheet(sv, name) is False:
            if DEBUG: print("   skip %s<-%s: both sheets fully read, neither prints the other's number" % (sv, name))
            continue
        for L, pairs in common_chain(vq, vp):
            strong = len(pairs) >= 3 or L >= 40.0
            R, t, rms, mx = rigid_fit(pairs)
            ov = 0.0 if ((sv in RAIL_PARCELS) != (name in RAIL_PARCELS)) else ov_area(unary_union([apply(p, R, t) for _, p in sheet]), anchor_geom)
            q_placed = unary_union([apply(p, R, t) for _, p in sheet])
            side_ok = True; side_note = ""; label_d = 0.0
            if name.upper() in NBPOS.get(sv, {}):          # where q's sheet writes p's number must land next to placed p
                lp = Point(apply(Point(NBPOS[sv][name.upper()]), R, t)); d = lp.distance(anchor_geom); side_ok &= d <= LABEL_REACH; side_note += " q-label->p %.0fm" % d; label_d = max(label_d, d)
            if sv.upper() in NBPOS.get(name, {}):          # where p's sheet writes q's number must land next to placed q
                lq = Point(apply(Point(NBPOS[name][sv.upper()]), Rp, tp)); d = lq.distance(q_placed); side_ok &= d <= LABEL_REACH; side_note += " p-label->q %.0fm" % d; label_d = max(label_d, d)
            r_ok, cross = rail_ok(sv, q_placed); side_note += " track-inside %.0fm%s" % (cross, "" if r_ok else " RAIL")
            theta_c = math.degrees(math.atan2(R[1, 0], R[0, 0])); dh = abs(((theta_c - prior_deg) + 180) % 360 - 180)
            heading_ok = dh <= HEADING_GATE_DEG or len(pairs) >= 3
            pc = puvi_centroid(sv); dpuvi = q_placed.centroid.distance(pc) if pc is not None else None; puvi_ok = True
            side_note += " heading %+.0f deg%s" % (dh, "" if heading_ok else " GATE") + ("" if dpuvi is None else " puvi-ref %.0fm" % dpuvi)
            ok = ((labelled and L >= 8.0) or strong) and ov <= max(2.0, 0.6 * L) and rms <= 1.5 and side_ok and heading_ok and puvi_ok and r_ok
            if DEBUG: print("   cand %s<-%s: L=%.1f pairs=%d labelled=%s rms=%.2f overlap=%.1f%s -> %s" % (sv, name, L, len(pairs), labelled, rms, ov, side_note, "OK" if ok else "rejected"))
            if not ok: continue
            cands.append({"name": name, "pairs": pairs, "L": L, "labelled": labelled, "R": R, "t": t, "rms": rms, "label_d": label_d,
                          "theta": math.degrees(math.atan2(R[1, 0], R[0, 0]))})
    if not cands: return None
    # anchors must agree: keep the largest group of candidates whose transforms agree within 1 m / 1 deg
    def agree(a, b):
        return abs(((a["theta"] - b["theta"]) + 180) % 360 - 180) <= 1.0 and float(np.linalg.norm(a["t"] - b["t"])) <= 2.0
    groups = []
    for c in cands:
        for g in groups:
            if all(agree(c, o) for o in g): g.append(c); break
        else: groups.append([c])
    sheet_area = unary_union([p for _, p in sheet]).area
    evals = []
    for group in groups:
        pairs_all = [pr for c in group for pr in c["pairs"]]
        R, t, rms, mx = rigid_fit(pairs_all); theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
        mine = unary_union([apply(p, R, t) for _, p in sheet])
        ov = sum(ov_area(mine, unary_union([apply(p, Rp, tp) for _, p in pp])) for nm, (pp, Rp, tp) in placed.items() if (sv in RAIL_PARCELS) == (nm in RAIL_PARCELS))
        tot_L = sum(c["L"] for c in group); limit = max(2.0, 0.6 * tot_L, 0.08 * sheet_area)   # older FMBs overlap later railway-land sheets; topology clips it
        lc = label_conflicts(sv, sheet, R, t, placed); r_ok, cross = rail_ok(sv, mine)
        base = sum(c["L"] + 10 * (len(c["pairs"]) - 2) + (15 if c["labelled"] else 0) - 0.5 * c["label_d"] for c in group) + 0.1 * len(group)
        score = base - 0.1 * ov - 20 * len(lc) + (0.5 * cross if sv in RAIL_PARCELS else -0.5 * cross)
        strong = any(len(c["pairs"]) >= 3 and c["labelled"] and c["label_d"] <= 15 for c in group)   # 3 vertex pairs + both numbers printed at the boundary
        valid = (ov <= limit or strong) and rms <= 1.5 and not lc and r_ok
        if DEBUG: print("   group %s: score %.1f rms=%.2f overlap with placed=%.1f m2 (limit %.1f) label conflicts=%s track-inside %.0fm -> %s" % ([c["name"] for c in group], score, rms, ov, limit, lc, cross, "valid" if valid else "invalid"))
        evals.append((score, valid, group, R, t, rms, mx, theta, mine, ov, tot_L))
    valid_evals = [e for e in evals if e[1]]
    if not valid_evals: return None                                 # every group lands on another sheet, breaks a label or straddles the track
    score, _, group, R, t, rms, mx, theta, mine, ov, tot_L = max(valid_evals, key=lambda e: e[0])
    f = lambda x, full, zero: max(0.0, min(1.0, (zero - x) / (zero - full)))
    chain_len = sum(c["L"] for c in group); anchors = [(c["name"], len(c["pairs"]) - 1, c["L"]) for c in group]
    conf = 40 * f(rms, 0.30, 1.5) + 30 * min(1.0, chain_len / 20.0) + 30 * min(2, len(group)) / 2
    if len(groups) > 1: conf -= 20                                # some candidate anchors disagreed
    if ov > max(2.0, 0.6 * tot_L, 0.08 * sheet_area): conf -= 10   # accepted on strong corroboration despite overlap
    supported = any(c["labelled"] for c in group)
    if not supported: conf -= 20                                  # no sheet prints the other's number: coincidental chains possible
    single_edge = len(group) == 1 and len(group[0]["pairs"]) == 2
    named = all(c["name"].upper() in mine_nb or sv.upper() in sheet_neighbours(c["name"]) for c in group)
    tentative_ok = single_edge and named
    if not supported: tentative_ok = False
    if any(c["name"] in TENTATIVE for c in group): conf = min(conf, 70)   # chained off a tentative anchor
    if single_edge: conf = min(conf, 79)                          # one edge only: never green
    return {"sheet": sheet, "R": R, "t": t, "theta": theta, "rms": rms, "max": mx, "anchors": anchors, "conf": round(max(0, conf)), "overlap": round(ov, 2),
            "rejected": [c["name"] for g in groups if g is not group for c in g], "single_edge": single_edge, "tentative_anchor": tentative_ok, "supported": supported}

def bearings_and_write(sv, sheet, R, t, method, conf, rms, anchors, out_dir):
    proj = Proj("EPSG:32644"); run = datetime.datetime.now().isoformat(timespec="seconds")
    rows, erows = [], []
    for props, g in sheet:
        pg = apply(g, R, t); assert abs(pg.area - g.area) < 1e-6 and abs(pg.length - g.length) < 1e-6
        p = dict(props); a = pg.area
        p.update({"area_sqm": round(a, 3), "area_are": round(a / 100, 4), "area_hect": round(a / 1e4, 6), "area_acre": round(a * ACRE, 5), "area_cent": round(a * ACRE * 100, 3),
                  "perimeter_m": round(pg.length, 3), "georef_method": method, "georef_confidence": conf, "georef_residual_m": round(rms, 3),
                  "anchored_to": ",".join(a_[0] for a_ in anchors), "georef_run": run})
        rows.append({**p, "geometry": pg})
        c = list(pg.exterior.coords); c = c[::-1] if Polygon(c).exterior.is_ccw else c       # clockwise
        start = max(range(len(c) - 1), key=lambda i: c[i][1]); c = c[start:-1] + c[:start]     # from northernmost vertex
        for i in range(len(c)):
            (x0, y0), (x1, y1) = c[i], c[(i + 1) % len(c)]; L = math.dist((x0, y0), (x1, y1))
            grid = (math.degrees(math.atan2(x1 - x0, y1 - y0))) % 360
            lon, lat = proj((x0 + x1) / 2, (y0 + y1) / 2, inverse=True); conv = proj.get_factors(lon, lat).meridian_convergence
            erows.append({"poly_id": props.get("poly_id"), "plot_no": props.get("plot_no"), "edge_no": i + 1, "length_m": round(L, 3), "bearing_grid_deg": round(grid, 2),
                          "bearing_true_deg": round((grid + conv) % 360, 2), "grid_convergence_deg": round(conv, 3), "geometry": LineString([(x0, y0), (x1, y1)])})
    out = os.path.join(out_dir, sv + "_parcels_AUTODEMO.gpkg")
    if os.path.exists(out): os.remove(out)
    gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:32644").to_file(out, layer="parcels", driver="GPKG")
    gpd.GeoDataFrame(erows, geometry="geometry", crs="EPSG:32644").to_file(out, layer="edges", driver="GPKG")
    return out

if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[1].startswith("35_"):
        set_village(sys.argv[1]); seed_sv = sys.argv[2]; targets = sys.argv[3].split(","); SEED_FILE = sys.argv[4] if len(sys.argv) > 4 else None
    else:
        seed_sv = sys.argv[1] if len(sys.argv) > 1 else "171"; targets = sys.argv[2].split(",") if len(sys.argv) > 2 else ["47A", "47B", "48A"]
    sheet, R, t, rms = seed_from_hand(seed_sv, SEED_FILE); prior = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    placed = {seed_sv: (sheet, R, t)}
    print("seed %s: rigid re-fit of the hand placement, rotation %.2f deg, hand-vs-rigid rms %.2f m" % (seed_sv, prior, rms))
    bearings_and_write(seed_sv, sheet, R, t, "seed", 100, rms, [], D)
    pending = list(targets); report = []
    while pending:
        best = None
        for sv in pending:
            r = place(sv, placed, prior)
            if r and (best is None or r["conf"] > best[1]["conf"]): best = (sv, r)
        if best is None: print("no more chains for", pending); break
        sv, r = best; pending.remove(sv)
        if not r.get("single_edge"): placed[sv] = (r["sheet"], r["R"], r["t"])
        elif r.get("tentative_anchor"):
            placed[sv] = (r["sheet"], r["R"], r["t"]); TENTATIVE.add(sv); print("   (%s placed on a single edge, named by its neighbour, heading agrees: TENTATIVE anchor)" % sv)
        else: print("   (%s placed on a single edge: written for review, NOT used as an anchor)" % sv)
        out = bearings_and_write(sv, r["sheet"], r["R"], r["t"], "chained", r["conf"], r["rms"], r["anchors"], D)
        # compare with the team's hand placement (pre-topology backup): centroid shift and mean vertex distance
        hb = glob.glob(os.path.join(BK, sv + "_parcels_modified*.gpkg")); mine = unary_union([apply(p, r["R"], r["t"]) for _, p in r["sheet"]])
        if hb: hand = unary_union(gpd.read_file(hb[0]).geometry.tolist()); shift = mine.centroid.distance(hand.centroid); haus = mine.hausdorff_distance(hand)
        else: shift = haus = float("nan")
        print("placed %-4s conf=%3d rms=%.2f m max=%.2f m rot=%.2f deg anchors=%s rejected=%s | vs team: centroid shift %.2f m, Hausdorff %.2f m | %s" % (sv, r["conf"], r["rms"], r["max"], r["theta"], [(a[0], a[1], round(a[2], 1)) for a in r["anchors"]], r["rejected"], shift, haus, os.path.basename(out)))
