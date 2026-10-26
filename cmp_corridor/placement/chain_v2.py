"""Chain placer v2 (south of 418): every sheet starts at its Puvi prior pose (prior.json: the Puvi survey polygon IS the sheet shape, rms 1-3 m,
so the prior gives position and orientation to about 10 m / 2 deg), then a rigid ICP closes it onto the already placed neighbours (outline points
within D of a neighbour outline are pulled onto it; D shrinks 25 -> 4 m), then the seams are closed by the smooth warp, T-junctions, despike and
the small polish, exactly as the chain did. Rigid sheets; nothing already placed moves. Usage: python chain_v2.py 420 423 187 ...
Writes <sv>_parcels_auto_V2.gpkg in this folder and a picture."""
import sys, math, json, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
from scipy.spatial import cKDTree
PRIOR = json.load(open("prior.json"))
NB = {
    "416": {"414": 60, "415": 60, "418": 250, "412": 20},
    "418": {"412": 20, "416": 250, "419": 170, "420": 126, "190": 100, "189": 40, "188": 43, "423": 19},
    "419": {"418": 170, "420": 90, "416": 10},
    "190": {"418": 100, "189": 60},
    "189": {"190": 60, "418": 40, "188": 22, "187": 62, "423": 10},
    "420": {"418": 126, "419": 90, "423": 308},
    "423": {"420": 308, "141": 86, "424": 78, "140": 75, "142": 53, "425": 32, "187": 31, "451": 26, "418": 19, "189": 10},
    "187": {"142": 85, "189": 62, "188": 59, "423": 31},
    "424": {"423": 78, "451": 41, "425": 27},
    "425": {"451": 135, "139": 100, "129": 68, "423": 32, "424": 27, "128": 19, "426": 14},
    "451": {"425": 135, "424": 41, "449": 30, "423": 26},
    "426": {"128": 118, "447": 103, "127": 88, "126": 72, "448": 72, "428": 50, "116": 41, "427": 40, "425": 14, "449": 9},
    "427": {"428": 141, "447": 62, "426": 40},
    "428": {"116": 166, "427": 141, "426": 50, "115": 22},
    "447": {"426": 103, "448": 73, "427": 62},
    "448": {"449": 110, "447": 73, "426": 72},
    "449": {"448": 110, "451": 30, "426": 9},
    "139": {"425": 100, "140": 72, "129": 68},
    "140": {"423": 75, "139": 72, "141": 31},
    "141": {"423": 86, "142": 35, "140": 31},
    "142": {"187": 85, "423": 53, "141": 35},
    "115": {"116": 149, "428": 22},
    "116": {"428": 166, "115": 149, "126": 77, "426": 41},
    "126": {"116": 77, "127": 74, "426": 72},
    "127": {"426": 88, "128": 84, "126": 74},
    "128": {"426": 118, "127": 84, "129": 47, "425": 19},
    "129": {"139": 68, "425": 68, "128": 47},
    "188": {"187": 59, "418": 43, "189": 22},
}
base = "D:/Projects/CMP_Roads_LPS/FMB_Vector/35_15_028"
ORDER = sys.argv[1:]
CLAMP = float(os.environ.get("CLAMP_DEG", "3.0"))  # a sheet may turn at most this much away from its prior orientation while closing onto its neighbours
# fixed base: the approved (green) sheets, his 414, and the option-A 415; everything red south of them is re-placed
BASE_FILES = {sv: f"{base}/{sv}_parcels_auto.gpkg" for sv in ("397", "398", "409", "410", "411", "412", "413", "415")}
BASE_FILES["414"] = f"{base}/414_parcels_modified.gpkg"
# sheets whose Puvi polygon is not the sheet (187: Puvi's 187 is a different 10 ha parcel): prior = the pose of the current delivered file
FILE_PRIOR = {"187"}  # 188 has a usable Puvi polygon (rms 2.1 m)
placed = {}
for sv, f in BASE_FILES.items():
    if sv in ORDER:
        continue
    g, U = load(f)
    placed[sv] = (g, U)
for f in sorted(os.listdir(".")):
    if f.endswith("_parcels_auto_V2.gpkg"):
        sv = f.split("_")[0]
        if sv not in placed and sv not in ORDER:
            g, U = load(f)
            placed[sv] = (g, U)


def dens_pts(geom, step=2.0):
    """densified exterior points of a polygon (or multipolygon) with their outward unit normals"""
    from shapely.geometry.polygon import orient
    pts, nrm = [], []
    for pp in getattr(geom, "geoms", [geom]):
        if pp.geom_type != "Polygon":
            continue
        L = orient(pp, 1.0).exterior  # counter-clockwise
        n = max(int(L.length / step), 8)
        for i in range(n):
            a = L.interpolate(i / n, normalized=True).coords[0]
            b = L.interpolate((i + 0.5) / n, normalized=True).coords[0]
            c = L.interpolate(((i - 0.5) % n) / n, normalized=True).coords[0]
            dx, dy = b[0] - c[0], b[1] - c[1]
            h = math.hypot(dx, dy) or 1.0
            pts.append(a)
            nrm.append((dy / h, -dx / h))  # outward for a CCW ring
    return np.array(pts), np.array(nrm)


def compose(Mi, M):
    """pose Mi applied after pose M (6-list affine)"""
    a, b, c, d, e, f = M
    a2, b2, c2, d2, e2, f2 = Mi
    return [a2 * a + b2 * c, a2 * b + b2 * d, c2 * a + d2 * c, c2 * b + d2 * d, a2 * e + b2 * f + e2, c2 * e + d2 * f + f2]


def rot_of(M):
    return math.degrees(math.atan2(M[2], M[0]))


def with_rotation(M, Us, deg):
    """same pose but with the sheet turned to `deg` about the centroid of its current image"""
    c = affine_transform(Us, M).centroid
    d = math.radians(deg - rot_of(M))
    ct, st = math.cos(d), math.sin(d)
    R = [ct, -st, st, ct, c.x - (ct * c.x - st * c.y), c.y - (st * c.x + ct * c.y)]
    return compose(R, M)


def icp(Us, M0, nb_polys, log, Ds=(25, 18, 12, 8, 5, 4, 4, 4)):
    """rigid refinement from pose M0: sheet outline points within D of a neighbour outline are matched to their nearest neighbour points,
    but only where the two edges face each other (outward normals opposed) - so a thin strip cannot be pulled onto a neighbour by its far
    side; Kabsch, repeat with a shrinking D; the rotation is clamped to +-CLAMP deg of the prior orientation (seams alone can leave a
    sheet free to pivot on a short contact)"""
    M = list(M0)
    th0 = rot_of(M0)
    B, NB_ = [], []
    for pp in nb_polys:
        b, nb = dens_pts(pp, 1.0)
        B.append(b)
        NB_.append(nb)
    B = np.vstack(B)
    NB_ = np.vstack(NB_)
    tB = cKDTree(B)
    hist = []

    def matches(M, D):
        """for every sheet outline point: the nearest neighbour-outline point whose edge faces it (the far edge of a thin strip is skipped)"""
        Up = affine_transform(Us, M)
        P, NP = dens_pts(Up, 1.5)
        K = min(40, len(B))
        d, j = tB.query(P, k=K)
        facing = np.einsum("ij,ikj->ik", NP, NB_[j]) < -0.3
        ok = (d < D) & facing
        first = np.argmax(ok, axis=1)
        sel = ok[np.arange(len(P)), first]
        return P, j[np.arange(len(P)), first], sel

    for D in Ds:
        P, j, sel = matches(M, D)
        if sel.sum() < 8:
            hist.append((D, int(sel.sum()), None, None))
            continue
        Mi, res, rot = rigid(P[sel], B[j[sel]])
        M = compose(Mi, M)
        dth = (rot_of(M) - th0 + 180) % 360 - 180
        if abs(dth) > CLAMP:
            M = with_rotation(M, Us, th0 + math.copysign(CLAMP, dth))
            P, j, sel = matches(M, D)  # translation-only re-fit at the clamped rotation
            if sel.sum() >= 8:
                t = (B[j[sel]] - P[sel]).mean(0)
                M = compose([1, 0, 0, 1, float(t[0]), float(t[1])], M)
                res = np.hypot(*(P[sel] + t - B[j[sel]]).T)
        hist.append((D, int(sel.sum()), round(float(np.sqrt((res ** 2).mean())), 2), round((rot_of(M) - th0 + 180) % 360 - 180, 2)))
    log.append(f"   icp (D, points, rms, turn from prior): {hist}")
    return M


ANCHOR = {"187": "189", "188": "187"}  # the neighbour the old chain hung these sheets on


def fit_to_file(p, g):
    """rigid fit of the PDF plots to a placed file (matched plot centroids by subdiv)"""
    g = g.copy()
    g["subdiv_no"] = g.subdiv_no.astype(str).str.strip()
    src, dst = [], []
    for _, r in p.iterrows():
        m = g[g.subdiv_no == str(r.subdiv_no).strip()]
        if str(r.subdiv_no).strip() not in ("", "-", "nan") and len(m) == 1:
            src.append(r.geometry.centroid.coords[0])
            dst.append(m.geometry.iloc[0].centroid.coords[0])
    M, res, rot = rigid(src, dst)
    return M, float(res.max())


def file_prior(sv, p):
    """prior for a sheet without a usable Puvi polygon: its old delivered pose, carried along with the motion of the neighbour it was
    hung on (old anchor file -> the anchor as placed now), so the old seam between them is kept as the starting point"""
    g_old = gpd.read_file(f"{base}/{sv}_parcels_auto.gpkg").to_crs(32644)
    M_old, res = fit_to_file(p, g_old)
    a = ANCHOR.get(sv)
    if a and a in placed:
        pa, _ = pdf_plots(a, verbose=False)
        Ma_old, _ = fit_to_file(pa, gpd.read_file(f"{base}/{a}_parcels_auto.gpkg").to_crs(32644))
        Ma_new, _ = fit_to_file(pa, placed[a][0])
        # motion old -> new of the anchor: Ma_new o Ma_old^-1
        a_, b_, c_, d_, e_, f_ = Ma_old
        det = a_ * d_ - b_ * c_
        inv = [d_ / det, -b_ / det, -c_ / det, a_ / det, (b_ * f_ - d_ * e_) / det, (c_ * e_ - a_ * f_) / det]
        Mt = compose(Ma_new, inv)
        return compose(Mt, M_old), res
    return M_old, res


def place(sv, log):
    p, sk = pdf_plots(sv, verbose=False)
    Us = make_valid(unary_union(list(p.geometry)))
    if sv in FILE_PRIOR:
        M0, fr = file_prior(sv, p)
        prior_note = f"prior = current file pose (fit residual {fr:.2f} m)"
    elif sv in PRIOR and PRIOR[sv]["rms_m"] < 15 and 0.6 < PRIOR[sv]["sheet_area"] / max(PRIOR[sv]["puvi_area"], 1) < 1.5:
        M0 = PRIOR[sv]["M"]
        prior_note = f"prior = Puvi polygon (rms {PRIOR[sv]['rms_m']} m)"
    else:
        log.append(f"{sv}: no usable prior")
        return None
    nbs = [o for o in NB.get(sv, {}) if o in placed]
    pool = nbs if nbs else [o for o in placed if placed[o][1].distance(affine_transform(Us, M0)) < 40]
    if not pool:
        log.append(f"{sv}: no placed neighbour near the prior pose; prior pose kept")
        M = M0
    else:
        nbP = [pp for o in pool for pp in getattr(placed[o][1], "geoms", [placed[o][1]]) if pp.geom_type == "Polygon"]
        M = icp(Us, M0, nbP, log)
    Up = affine_transform(Us, M)
    c0, c1 = affine_transform(Us, M0).centroid, Up.centroid
    MAXMOVE = float(os.environ.get("MAXMOVE", "25"))
    if math.hypot(c1.x - c0.x, c1.y - c0.y) > MAXMOVE:
        # a rigid fit on too little contact can run away; keep the prior pose and let the seam closing do the rest
        log.append(f"   icp moved the sheet {math.hypot(c1.x - c0.x, c1.y - c0.y):.1f} m from its prior (> {MAXMOVE} m): rejected, prior pose kept")
        M = list(M0)
        Up = affine_transform(Us, M)
        c1 = Up.centroid
    rot0 = math.degrees(math.atan2(M0[2], M0[0]))
    rot1 = math.degrees(math.atan2(M[2], M[0]))
    log.append(f"{sv}: {prior_note}; moved from the prior by {math.hypot(c1.x - c0.x, c1.y - c0.y):.1f} m, {((rot1 - rot0 + 180) % 360 - 180):+.2f} deg; neighbours used {pool}")
    ov = {o: round(Up.intersection(placed[o][1]).area, 1) for o in placed if Up.intersection(placed[o][1]).area > 1}
    if sum(ov.values()) > 0.05 * Us.area:
        log.append(f"{sv}: STOP - overlaps {ov} m2 after the rigid fit")
        return None
    g = p.copy()
    g["geometry"] = [affine_transform(x, M) for x in p.geometry]
    g = g.set_crs(32644, allow_override=True)
    paths = []
    Ug = make_valid(unary_union(list(g.geometry)))
    for o in placed:
        if Ug.distance(placed[o][1]) < 3.0 or Ug.intersects(placed[o][1]):
            Ro = ring(placed[o][1])
            paths.append(Ro + [Ro[0]])
    mv = []
    if paths:
        g, mv = warp_snap(g, paths, tol=3.0)
    dlog0 = []
    g["geometry"] = [x if x.is_valid else despike(x, dlog0) for x in g.geometry]
    if dlog0:
        log.append(f"   despiked after the warp: {dlog0}")
    Ug = make_valid(unary_union(list(g.geometry)))
    # option A squeeze: where the rigid sheet still contradicts a placed neighbour (overlap beyond the 3 m seam closing), its outline is
    # warped further onto that neighbour (up to 8 m); whatever is left is clipped - the placed side never moves
    big = {o: Ug.intersection(placed[o][1]).area for o in placed if Ug.intersection(placed[o][1]).area > 25.0}
    if big:
        cut = unary_union([placed[o][1] for o in big])
        newg, lost = [], []
        for k, x in zip(g.subdiv_no.astype(str), g.geometry):
            if x.intersection(cut).area > 0.01:
                y = as_poly(make_valid(x.difference(cut)))
                if y.is_empty or y.area < 1.0:
                    lost.append((k, round(x.area, 1)))
                    y = x
                newg.append(y)
            else:
                newg.append(x)
        g["geometry"] = newg
        log.append(f"   squeeze: {sum(big.values()):.1f} m2 overlapping {list(big)} clipped from {sv} (the placed side never moves)" + (f"; plots entirely inside the neighbour, kept unclipped: {lost}" if lost else ""))
        Ug = make_valid(unary_union(list(g.geometry)))
    shared = {o: round(Ug.boundary.intersection(placed[o][1].boundary.buffer(0.05)).length) for o in placed if Ug.distance(placed[o][1]) < 0.1}
    log.append(f"   seams: {len(mv)} outline vertices moved, max {max(mv) if mv else 0} m; shared boundary now {shared} m (Puvi expects {NB.get(sv, {})})")
    g = tjunctions(g, [placed[o][0] for o in placed if Ug.distance(placed[o][1]) < 3.0])
    dlog = []
    g["geometry"] = [x if x.is_valid else despike(x, dlog) for x in g.geometry]
    if dlog:
        log.append(f"   despiked: {dlog}")
    Ug = make_valid(unary_union(list(g.geometry)))
    for o in placed:
        inter = Ug.intersection(placed[o][1])
        if 0.01 < inter.area < 25.0:
            g["geometry"] = [make_valid(x.difference(placed[o][1])) if x.intersection(placed[o][1]).area > 0.005 else x for x in g.geometry]
            log.append(f"   polish: {inter.area:.2f} m2 overlap with {o} cut from {sv}")
    Ug = make_valid(unary_union(list(g.geometry)))
    Uall = make_valid(unary_union([Ug] + [placed[o][1] for o in placed if placed[o][1].distance(Ug) < 1]))
    for pp in getattr(Uall, "geoms", [Uall]):
        for r_ in pp.interiors:
            h = Polygon(r_)
            if 0.01 < h.area < 25.0 and h.distance(Ug) < 0.001:
                j = min(range(len(g)), key=lambda i: g.geometry.iloc[i].distance(h))
                g.loc[g.index[j], "geometry"] = make_valid(unary_union([g.geometry.iloc[j], h]))
                log.append(f"   polish: {h.area:.2f} m2 gap filled into {sv}/{g.subdiv_no.iloc[j]}")
    g["geometry"] = [x if x.geom_type == "Polygon" else max(x.geoms, key=lambda y: y.area) for x in g.geometry]
    # small plots can collapse in the seam closing (all their vertices snap onto one path): restore the rigid shape, minus what is placed
    rigid_geoms = [affine_transform(x, M) for x in p.geometry]
    placed_union = unary_union([placed[o][1] for o in placed if placed[o][1].distance(Ug) < 5])
    fixed = []
    for i, (x, r0) in enumerate(zip(g.geometry, rigid_geoms)):
        if x is None or x.is_empty or x.area < 0.5 * r0.area:
            y = as_poly(make_valid(r0.difference(placed_union)))
            others = unary_union([z for j, z in enumerate(g.geometry) if j != i and z is not None and not z.is_empty])
            y = as_poly(make_valid(y.difference(others)))
            g.loc[g.index[i], "geometry"] = y
            fixed.append((str(g.subdiv_no.iloc[i]), round(r0.area, 1), round(0 if x is None or x.is_empty else x.area, 1), round(y.area, 1)))
    if fixed:
        log.append(f"   small plots restored to their rigid shape (plot, sheet m2, after seams m2, restored m2): {fixed}")
    Ug = make_valid(unary_union(list(g.geometry)))
    ovs = {o: round(Ug.intersection(placed[o][1]).area, 2) for o in placed if Ug.intersection(placed[o][1]).area > 0.01}
    bad = [k for k, x in zip(kd(g), g.geometry) if not x.is_valid]
    sa = {rr.subdiv_no: rr.geometry.area for _, rr in sk.iterrows() if rr.subdiv_no}
    dev = sorted([(rr.subdiv_no, round((rr.geometry.area / sa[rr.subdiv_no] - 1) * 100, 1)) for _, rr in g.iterrows() if rr.subdiv_no in sa], key=lambda t: abs(t[1]))
    log.append(f"   after seams: overlaps {ovs}, invalid {bad}, plot changes (largest) {dev[-3:]}")
    ok = not ovs and not bad
    o = finish(g, sv, f"auto 2026-10-22 v2: rigid PDF sheet started at its Puvi position, closed onto the placed neighbours {pool} by a rigid fit (turn clamped to {CLAMP} deg), seams closed by the smooth warp (max {max(mv) if mv else 0} m)")
    return o, ok


if __name__ == "__main__":
    log = []
    done = []
    for sv in ORDER:
        try:
            res = place(sv, log)
        except Exception as e:
            import traceback
            print("\n".join(log))
            print(f"CRASH at {sv}: {e}")
            traceback.print_exc()
            raise SystemExit(1)
        if res is None:
            print("\n".join(log))
            print(f"STOP at {sv}")
            break
        o, ok = res
        o.to_crs(4326).to_file(f"{sv}_parcels_auto_V2.gpkg", layer=f"{sv}_parcels_auto", driver="GPKG")
        # the 4326 round trip can turn a near-degenerate spike invalid: read back, repair, rewrite
        rb = gpd.read_file(f"{sv}_parcels_auto_V2.gpkg").to_crs(32644)
        if not rb.geometry.is_valid.all():
            fl = []
            rb["geometry"] = [x if x.is_valid else despike(x, fl) for x in rb.geometry]
            rb.to_crs(4326).to_file(f"{sv}_parcels_auto_V2.gpkg", layer=f"{sv}_parcels_auto", driver="GPKG")
            log.append(f"   repaired after the write round trip: {fl}")
            o = rb
        placed[sv] = (o, make_valid(unary_union(list(o.geometry))))
        done.append(sv)
        log.append(f"   written {sv}_parcels_auto_V2.gpkg ({len(o)} plots, {round(o.area_acre.sum(), 3)} ac) {'OK' if ok else 'CHECK'}")
        if not ok:
            print("\n".join(log))
            print(f"STOP after {sv}: needs a look")
            break
    else:
        print("\n".join(log))
    seams({sv: placed[sv][0] for sv in placed}, "ALL placed")
    if done:
        picture("chain_v2_picture.png", "v2: orange = base | red = re-placed " + ", ".join(done), [(placed[o][0], "ctx") for o in placed if o not in done] + [(placed[o][0], (255, 60, 60)) for o in done], sat_file="sat_big.png", X0=409650.0, Y1=1421800.0, pad=40)
