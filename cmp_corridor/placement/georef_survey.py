"""Rigid sheet for any corridor survey: PDF line work -> plots (labels from the sketch vector), rigid pose (rotation + shift, true sheet scale)
fitted to the CMP fine-tuned parcels of that survey (their position is roughly right, their shape is the sheet conformed to Puvi), written as
<vc>_<sv>_rigid.gpkg in this folder. Prints the fit quality. Usage: python georef_survey.py <village_code> <survey> [<survey> ...]"""
import sys, os, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
from scipy.spatial import cKDTree
from scipy.optimize import minimize
VC = sys.argv[1]
MERGED = "D:/Projects/CMP_Roads_LPS/FMB_Georeferenced/CMP_my_Part.geojson"
M = gpd.read_file(MERGED).to_crs(32644)
M["kide"] = M.kide.astype(str)
M["village_code"] = M.village_code.astype(str)
M["survey"] = M.kide.str.split("/").str[0]


def sheet_plots(vc, sv):
    pg = fitz.open(f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{vc}/{sv}.pdf")[0]
    W, H = pg.rect.width, pg.rect.height
    lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2)
        col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1":
            continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    polys = [p for p in polygonize(unary_union(lines)) if p.area > 1.0]
    sk = gpd.read_file(SK, layer=f"{vc}_{sv}")
    sk["subdiv_no"] = sk.subdiv_no.astype(str).str.strip()
    big, refp = max(polys, key=lambda p: p.area), max(sk.geometry, key=lambda g: g.area)
    polys = [translate(p, refp.centroid.x - big.centroid.x, refp.centroid.y - big.centroid.y) for p in polys]
    rows = []
    for p in polys:
        hit = [r.subdiv_no for _, r in sk.iterrows() if r.geometry.intersection(p).area > 0.6 * p.area and r.subdiv_no]
        rows.append({"subdiv_no": hit[0] if hit else "-", "geometry": p})
    return gpd.GeoDataFrame(rows, geometry="geometry"), sk


def dens(g, step=3.0):
    L = (g if g.geom_type == "Polygon" else max(g.geoms, key=lambda x: x.area)).exterior
    n = max(int(L.length / step), 8)
    return np.array([L.interpolate(i / n, normalized=True).coords[0] for i in range(n)])


def fit(Us, Uc, pairs=None):
    """rigid pose of the sheet union onto the CMP union: symmetric outline distance, 12 starts; if matched plot centroids exist they are added"""
    A = dens(Us); cA = np.array(Us.centroid.coords[0]); A0 = A - cA
    B = dens(Uc); cB = np.array(Uc.centroid.coords[0]); tB = cKDTree(B)
    P0 = np.array([p[0] for p in pairs]) - cA if pairs else None
    Q0 = np.array([p[1] for p in pairs]) if pairs else None

    def cost(x):
        th, tx, ty = x
        ct, st = math.cos(th), math.sin(th)
        R = np.array([[ct, st], [-st, ct]])
        P = A0 @ R + np.array([tx, ty])
        c = (np.sqrt((tB.query(P)[0] ** 2).mean()) + np.sqrt((cKDTree(P).query(B)[0] ** 2).mean())) / 2
        if P0 is not None and len(P0) >= 2:
            c += np.sqrt((((P0 @ R + np.array([tx, ty])) - Q0) ** 2).sum(1).mean())
        return c
    best = None
    for th0 in np.linspace(0, 2 * math.pi, 12, endpoint=False):
        r = minimize(cost, [th0, cB[0], cB[1]], method="Nelder-Mead", options={"xatol": 0.02, "fatol": 0.01, "maxiter": 1500, "initial_simplex": [[th0, cB[0], cB[1]], [th0 + 0.15, cB[0], cB[1]], [th0, cB[0] + 15, cB[1]], [th0, cB[0], cB[1] + 15]]})
        if best is None or r.fun < best.fun:
            best = r
    th, tx, ty = best.x
    ct, st = math.cos(th), math.sin(th)
    return [ct, -st, st, ct, tx - (ct * cA[0] - st * cA[1]), ty - (st * cA[0] + ct * cA[1])], best.fun, math.degrees(th) % 360


for sv in sys.argv[2:]:
    try:
        p, sk = sheet_plots(VC, sv)
    except Exception as e:
        print(f"{VC} {sv}: sheet not usable ({e})")
        continue
    Us = make_valid(unary_union(list(p.geometry)))
    c = M[(M.village_code == VC) & (M.survey == sv) & (~M.source.astype(str).str.startswith("georef"))].copy()
    if not len(c):
        print(f"{VC} {sv}: no CMP parcels to start from")
        continue
    c["subdiv_no"] = c.subdiv_no.astype(str).str.strip()
    Uc = make_valid(unary_union(list(c.geometry)))
    pairs = []
    for s in set(p.subdiv_no) & set(c.subdiv_no) - {"-", ""}:
        a = p[p.subdiv_no == s]; b = c[c.subdiv_no == s]
        if len(a) == 1 and len(b) == 1:
            pairs.append((a.geometry.iloc[0].centroid.coords[0], b.geometry.iloc[0].centroid.coords[0]))
    Mx, rms, rot = fit(Us, Uc, pairs if len(pairs) >= 3 else None)
    g = p.copy()
    g["geometry"] = [affine_transform(x, Mx) for x in p.geometry]
    g = g.set_crs(32644, allow_override=True)
    g["kide"] = [f"{sv}/{s}" if s not in ("-", "") else sv for s in g.subdiv_no]
    g["survey_no"] = sv; g["village_code"] = VC; g["poly_id"] = range(1, len(g) + 1)
    g["area_sqm"] = g.geometry.area.round(1); g["area_acre"] = (g.geometry.area / 4046.856).round(4)
    out = f"{VC}_{sv}_rigid.gpkg"
    g.to_crs(4326).to_file(out, layer=f"{sv}_parcels_auto", driver="GPKG")
    Up = make_valid(unary_union(list(g.geometry)))
    print(f"{VC} {sv}: sheet {len(p)} plots {Us.area:.0f} m2 (CMP {len(c)} plots {Uc.area:.0f} m2, ratio {Uc.area / Us.area:.2f}) | pose fit rms {rms:.1f} m on {len(pairs)} matched plots, centroid shift {Up.centroid.distance(Uc.centroid):.1f} m | unnamed {(p.subdiv_no == '-').sum()} -> {out}")
