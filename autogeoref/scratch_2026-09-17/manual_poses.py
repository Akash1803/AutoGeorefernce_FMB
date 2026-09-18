"""Recover a rigid pose (rotation + shift, scale 1) for every survey in a hand-georeferenced file,
robust to vertex edits: align polygon centroids by poly_id first, then refine on nearest vertices."""
import sys, os, json, math, pickle, numpy as np, geopandas as gpd
from shapely.geometry import shape
from shapely.ops import unary_union
from shapely import affinity
VILLAGE = sys.argv[1]; MANUAL = sys.argv[2]; D = "D:/Projects/Tambaram_Chengalpattu/FMB_Vector/" + VILLAGE
man = gpd.read_file(MANUAL)
if man.crs and man.crs.to_epsg() != 32644: man = man.to_crs(32644)
sk = man["sketch_id"].astype(str).str.strip().str.replace(VILLAGE + "_", "", regex=False)
man["sv"] = np.where(man["sketch_id"].notna(), sk, np.where(man["survey_no"].notna(), man["survey_no"].astype(str).str.strip(), "?"))
def rigid(P, Q):
    P, Q = np.asarray(P, float), np.asarray(Q, float); pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc); U, S, Vt = np.linalg.svd(H); R = Vt.T @ U.T
    if np.linalg.det(R) < 0: Vt[-1] *= -1; R = Vt.T @ U.T
    t = qc - R @ pc; res = np.linalg.norm((R @ P.T).T + t - Q, axis=1); scale = S.sum() / ((P - pc) ** 2).sum()
    return R, t, float(np.sqrt((res ** 2).mean())), float(res.max()), float(scale)
poses = {}
for sv, grp in man[man["sv"] != "?"].groupby("sv"):
    p = os.path.join(D, sv + "_parcels.geojson")
    if not os.path.exists(p): print("  %-4s no sheet file" % sv); continue
    sheet = {int(f["properties"]["poly_id"]): shape(f["geometry"]) for f in json.load(open(p, encoding="utf-8"))["features"]}
    mg = {}
    for r in grp.itertuples():
        try: mg[int(r.poly_id)] = r.geometry
        except Exception: pass
    common = [k for k in sheet if k in mg]
    if len(common) < 2:
        # single polygon: use its vertices via nearest matching after centroid alignment (rotation search)
        common = list(sheet)
    # stage 1: rigid fit on polygon centroids (poly_id correspondence), or rotation search for 1 polygon
    if len(common) >= 2:
        R, t, _, _, _ = rigid([sheet[k].centroid.coords[0] for k in common], [mg[k].centroid.coords[0] for k in common])
    else:
        k = common[0]; best = None
        for deg in np.arange(0, 360, 0.5):
            gs = affinity.rotate(sheet[k], deg, origin="centroid"); gs = affinity.translate(gs, mg[k].centroid.x - gs.centroid.x, mg[k].centroid.y - gs.centroid.y)
            d = gs.symmetric_difference(mg[k]).area
            if best is None or d < best[0]: best = (d, deg)
        th = math.radians(best[1]); R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
        c = np.array(sheet[k].centroid.coords[0]); t = np.array(mg[k].centroid.coords[0]) - R @ c
    # stage 2: refine on nearest vertices (sheet vertex -> nearest manual vertex of the same polygon, within 4 m)
    rms = mx = scale = float('nan'); P = []
    for _ in range(3):
        P, Q = [], []
        for k in common:
            mv = np.array(list((mg[k].geoms[0] if mg[k].geom_type == "MultiPolygon" else mg[k]).exterior.coords)[:-1])
            for v in list(sheet[k].exterior.coords)[:-1]:
                w = R @ np.array(v) + t; d = np.linalg.norm(mv - w, axis=1); j = int(d.argmin())
                if d[j] <= 4.0: P.append(v); Q.append(mv[j])
        if len(P) < 3: break
        R, t, rms, mx, scale = rigid(P, Q)
    theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    # implied stretch: best similarity scale on the same pairs
    poses[sv] = {"R": R, "t": t, "theta": theta, "rms": rms, "max": mx, "scale": scale, "n_pairs": len(P), "n_polys": len(sheet)}
    print("  %-4s polys %3d | heading %7.2f | implied scale %.4f | rigid rms %.2f m (max %.2f) on %d vertex pairs" % (sv, len(sheet), theta, scale, rms, mx, len(P)))
pickle.dump(poses, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "manual_poses_%s.pkl" % VILLAGE), "wb"))
thetas = [v["theta"] for v in poses.values()]; print("headings: median %.2f, spread %s" % (float(np.median(thetas)), [round(x, 1) for x in sorted(thetas)]))
