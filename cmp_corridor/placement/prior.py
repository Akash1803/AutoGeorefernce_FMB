"""Puvi prior pose for every corridor sheet of Ponmar: rigid fit (rotation + shift, true sheet scale) of the PDF sheet outline to the Puvi survey
polygon, by symmetric boundary distance from several starts. Then every delivered sheet is compared with its prior (rotation and centroid
difference) so the drift of the chain can be measured against an independent frame. Writes prior.json and prior_check.csv."""
import sys, os, math, json, glob
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
from scipy.optimize import minimize
import pandas as pd
PUVI = "D:/Data/Puvi_data/vector/35/15/028/vector/35_15_028_chengalpattu_vandalur_ponmar_vector.shp"
base = "D:/Projects/CMP_Roads_LPS/FMB_Vector/35_15_028"
pv = gpd.read_file(PUVI).to_crs(32644)
pv["survey_no"] = pv.survey_no.astype(str).str.strip()


def dens(geom, step=4.0):
    L = geom.exterior if geom.geom_type == "Polygon" else max(geom.geoms, key=lambda g: g.area).exterior
    n = max(int(L.length / step), 8)
    return np.array([L.interpolate(i / n, normalized=True).coords[0] for i in range(n)])


def pose_M(th, tx, ty, c):
    """rigid pose: rotate the sheet (about its centroid c) by th radians, then shift so the centroid lands on (tx, ty)"""
    ct, st = math.cos(th), math.sin(th)
    return [ct, -st, st, ct, tx - (ct * c[0] - st * c[1]), ty - (st * c[0] + ct * c[1])]


def fit_prior(Us, Pv):
    """rigid pose of the sheet outline on the Puvi polygon: symmetric nearest-point RMS between densified boundaries (KD-trees), 12 starts"""
    from scipy.spatial import cKDTree
    A = dens(Us, 3.0); c = np.array(Us.centroid.coords[0]); A0 = A - c
    Bp = dens(Pv, 3.0); cp = np.array(Pv.centroid.coords[0])
    tB = cKDTree(Bp)

    def cost(x):
        th, tx, ty = x
        ct, st = math.cos(th), math.sin(th)
        P = A0 @ np.array([[ct, st], [-st, ct]]) + np.array([tx, ty])
        d1 = tB.query(P)[0]
        d2 = cKDTree(P).query(Bp)[0]
        return (np.sqrt((d1 ** 2).mean()) + np.sqrt((d2 ** 2).mean())) / 2

    def axis(g):
        r = g.minimum_rotated_rectangle
        xs = list(r.exterior.coords)
        e = max(((math.dist(xs[i], xs[i + 1]), i) for i in range(4)))[1]
        return math.atan2(xs[e + 1][1] - xs[e][1], xs[e + 1][0] - xs[e][0])
    da = axis(Pv) - axis(Us)
    starts = [da + k * math.pi / 2 for k in range(4)] + [k * math.pi / 4 for k in range(8)]
    best = None
    for th0 in starts:
        r = minimize(cost, [th0, cp[0], cp[1]], method="Nelder-Mead", options={"xatol": 0.02, "fatol": 0.01, "maxiter": 1500, "initial_simplex": [[th0, cp[0], cp[1]], [th0 + 0.15, cp[0], cp[1]], [th0, cp[0] + 15, cp[1]], [th0, cp[0], cp[1] + 15]]})
        if best is None or r.fun < best.fun:
            best = r
    th, tx, ty = best.x
    return {"theta_deg": round(math.degrees(th) % 360, 2), "cx": round(tx, 2), "cy": round(ty, 2), "rms_m": round(float(best.fun), 2), "M": pose_M(th, tx, ty, c)}


def current_pose(sv, p, g):
    """rotation and centroid of the delivered sheet relative to the sheet frame: rigid fit of matched plot centroids (by subdiv), else outline vertices"""
    p = p.copy(); g = g.copy()
    p["subdiv_no"] = p.subdiv_no.astype(str).str.strip(); g["subdiv_no"] = g.subdiv_no.astype(str).str.strip()
    src, dst = [], []
    for _, r in p.iterrows():
        m = g[g.subdiv_no == r.subdiv_no]
        if r.subdiv_no not in ("", "-", "nan") and len(m) == 1:
            src.append(r.geometry.centroid.coords[0]); dst.append(m.geometry.iloc[0].centroid.coords[0])
    if len(src) < 2:
        Us = make_valid(unary_union(list(p.geometry))); Ug = make_valid(unary_union(list(g.geometry)))
        rs, rg = ring(Us), ring(Ug)
        # match by the longest edge
        def longest(r):
            n = len(r); k = max(range(n), key=lambda i: math.dist(r[i], r[(i + 1) % n])); return [r[k], r[(k + 1) % n]]
        src, dst = longest(rs), longest(rg)
    M, res, rot = rigid(src, dst)
    Us = make_valid(unary_union(list(p.geometry)))
    cg = affine_transform(Us, M).centroid
    return rot, cg.x, cg.y, float(res.max())


if __name__ == "__main__":
    svs = sys.argv[1:] or ["397", "398", "409", "410", "411", "412", "413", "414", "415", "416", "418", "419", "420", "190", "189", "187", "188", "423", "142", "141", "140", "424", "425", "451", "139", "129", "128", "426", "127", "126", "116", "115", "428", "427", "447", "448", "449"]
    out = json.load(open("prior.json")) if os.path.exists("prior.json") else {}
    rows = []
    for sv in svs:
        P = pv[pv.survey_no == sv]
        if not len(P):
            print(sv, "not in Puvi"); continue
        Pv = make_valid(unary_union(list(P.geometry)))
        p, sk = pdf_plots(sv, verbose=False)
        Us = make_valid(unary_union(list(p.geometry)))
        pr = fit_prior(Us, Pv)
        pr["puvi_area"] = round(Pv.area); pr["sheet_area"] = round(Us.area)
        out[sv] = pr
        row = {"survey": sv, "prior_theta": pr["theta_deg"], "prior_rms_m": pr["rms_m"], "area_ratio": round(Us.area / Pv.area, 2)}
        f = f"{base}/{sv}_parcels_auto.gpkg"
        if not os.path.exists(f): f = f"{base}/{sv}_parcels_modified.gpkg"
        if os.path.exists(f):
            g = gpd.read_file(f).to_crs(32644)
            rot, cx, cy, res = current_pose(sv, p, g)
            dth = (rot - pr["theta_deg"] + 180) % 360 - 180
            row.update({"cur_theta": round(rot, 2), "d_theta": round(dth, 2), "d_pos_m": round(math.hypot(cx - pr["cx"], cy - pr["cy"]), 1), "dx": round(cx - pr["cx"], 1), "dy": round(cy - pr["cy"], 1), "fit_res": round(res, 2)})
        rows.append(row)
        print(row, flush=True)
    json.dump(out, open("prior.json", "w"), indent=1)
    df = pd.DataFrame(rows); pd.set_option("display.width", 200)
    print(df.to_string(index=False))
    df.to_csv("prior_check.csv", index=False)
