"""Annotation check, read 1 (machine): for every survey of CMP_my_Part with a PDF sheet
 - PDF plots = polygonized line work (sheet frame X = px/PT, Y = -py/PT)
 - printed plot numbers = read5 (glyph outlines, plot-number size separated from dimensions) -> in the PDF plot that holds them,
   else the nearest plot within 6 m (numbers printed beside thin strips)
 - layer plots of the survey <-> PDF plots: rigid fit of the sheet outline onto the layer outline (shape only, no labels used),
   then one-to-one by overlap (IoU)
 - verdict per layer plot: ok | WRONG (PDF prints another number) | MISSING-NUMBER (layer has none, PDF prints one)
   | NOT-PRINTED (layer has a number, PDF plot prints none) | UNMATCHED (no PDF plot overlaps well)
-> qc/ann_read1.csv (one row per layer plot), qc/ann_pdfplots.csv (one row per PDF plot), qc/ann_sheets.csv"""
import sys, os, math, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "qc", "r5")); sys.path.insert(1, r"D:\code\FMB_to_GeoJSON"); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, geopandas as gpd, fitz
import read5
from shapely import make_valid
from shapely.geometry import LineString, Point
from shapely.ops import unary_union, polygonize
from shapely.affinity import affine_transform
from scipy.spatial import cKDTree
from scipy.optimize import minimize, linear_sum_assignment
PT = 72 / 25.4 * 2          # pt per metre at 1:500
SK = "D:/Projects/CMP_Roads_LPS/FMB_Sketches"
L = gpd.read_file(os.environ.get("ANN_FILE", r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_my_Part.geojson")).to_crs(32644)
PFX = os.environ.get("ANN_PFX", "ann")
L["kide"] = L.kide.astype(str); L["vc"] = L.village_code.astype(str); L["sv"] = L.kide.str.split("/").str[0]
L["lab"] = L.subdiv_no.astype(str).str.strip().replace({"None": "", "nan": "", "-": ""})
only = set(sys.argv[1:])
T = pd.read_excel("D:/Projects/CMP_Roads_LPS/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
REG = T.groupby([T.vc, T.fmb_survey_no.str.upper().str.replace(" ", "")]).sub_division.apply(lambda s: {"".join(ch for ch in x.upper() if ch.isalnum()) for x in s} - {""}).to_dict()
norm = lambda s: "".join(ch for ch in str(s).upper() if ch.isalnum())


def pdf_polys(pg):
    W, H = pg.rect.width, pg.rect.height; lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2); col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1": continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    return [p for p in polygonize(unary_union(lines)) if p.area > 0.5]


def dens(g, step):
    g = g if g.geom_type == "Polygon" else max(g.geoms, key=lambda q: q.area)
    n = max(int(g.exterior.length / step), 16)
    return np.array([g.exterior.interpolate(i / n, normalized=True).coords[0] for i in range(n)])


def fit(Us, Uc):
    step = max(1.0, math.sqrt(Uc.area) / 60)
    A = dens(Us, step); cA = np.array(Us.centroid.coords[0]); A0 = A - cA; B = dens(Uc, step); tB = cKDTree(B); cB = np.array(Uc.centroid.coords[0])
    def cost(x):
        th, tx, ty = x; ct, st = math.cos(th), math.sin(th); P = A0 @ np.array([[ct, st], [-st, ct]]) + [tx, ty]
        return (np.sqrt((tB.query(P)[0] ** 2).mean()) + np.sqrt((cKDTree(P).query(B)[0] ** 2).mean())) / 2
    best = min((minimize(cost, [t, *cB], method="Nelder-Mead", options={"maxiter": 3000, "xatol": 1e-4, "fatol": 1e-4}) for t in np.linspace(0, 2 * math.pi, 16, endpoint=False)), key=lambda r: r.fun)
    th, tx, ty = best.x; ct, st = math.cos(th), math.sin(th)
    return [ct, -st, st, ct, tx - (ct * cA[0] - st * cA[1]), ty - (st * cA[0] + ct * cA[1])], best.fun


rows, prow, srow = [], [], []; WARPED = []; PAIRS = {}; RIGID = {}
for (vc, sv), cur in L.groupby(["vc", "sv"]):
    if only and f"{vc}:{sv}" not in only: continue
    pdf = f"{SK}/{vc}/{sv}.pdf"
    if not os.path.exists(pdf):
        srow.append({"village": vc, "survey": sv, "status": "no PDF", "layer_plots": len(cur)})
        for r in cur.itertuples(): rows.append({"village": vc, "survey": sv, "kide": r.kide, "layer_label": r.lab, "verdict": "NO PDF"})
        continue
    doc = fitz.open(pdf); pg = doc[0]; H = pg.rect.height
    polys = pdf_polys(pg)
    try:
        labs, info = read5.read(pdf, 500); ok = {id(a) for a in read5.plot_numbers(labs)}
    except Exception as e:
        labs, info, ok = [], {}, set()
    nums = [(a["text"], Point(a["px"] / PT, -a["py"] / PT)) for a in labs if id(a) in ok]
    # numbers -> PDF plots
    reg = REG.get((vc, sv.upper()), set())
    if reg: nums = [(t, p) for t, p in nums if norm(t) in reg]          # only numbers the A-register knows for this survey
    pn = {k: [] for k in range(len(polys))}
    for t, p in nums:
        for k, q in enumerate(polys):
            d = q.distance(p)
            if d == 0: pn[k].append(t)
            elif d <= 3.0: pn[k].append(t + "*")                            # printed beside (thin plots)
    if not polys:
        srow.append({"village": vc, "survey": sv, "status": "PDF line work does not close", "layer_plots": len(cur), "numbers": len(nums)})
        for r in cur.itertuples(): rows.append({"village": vc, "survey": sv, "kide": r.kide, "layer_label": r.lab, "verdict": "NO PDF PLOTS"})
        continue
    Us = make_valid(unary_union(polys)); Uc = make_valid(unary_union(list(cur.geometry)))
    if len(cur) == 1 and len(polys) == 1:
        Mx, rms = None, 0.0; W_ = [polys[0]]
    else:
        Mx, rms = fit(Us, Uc); W_ = [affine_transform(q, Mx) for q in polys]; RIGID[(vc, sv)] = Mx
        # distorted team vectors: anchor on unique label pairs and warp the sheet locally (thin-plate spline), outliers dropped
        lab_sheet = {}
        for k in range(len(polys)):
            ins = sorted({norm(t) for t in pn[k] if not t.endswith("*")})
            if len(ins) == 1: lab_sheet.setdefault(ins[0], []).append(k)
        lab_lay = {}
        for i, r in enumerate(cur.itertuples()): lab_lay.setdefault(norm(r.lab), []).append(i)
        pr = [(polys[ks[0]].representative_point().coords[0], cur.geometry.iloc[ls[0]].representative_point().coords[0]) for t, ks in lab_sheet.items() if t and len(ks) == 1 and len(lab_lay.get(t, [])) == 1 for ls in [lab_lay[t]]]
        if len(pr) >= 6:
            from scipy.interpolate import RBFInterpolator
            Pq = np.array([p[0] for p in pr]); Qq = np.array([p[1] for p in pr]); keep = np.ones(len(pr), bool)
            for _ in range(4):
                tps = RBFInterpolator(Pq[keep], Qq[keep], kernel="thin_plate_spline", smoothing=2.0)
                res = np.hypot(*(tps(Pq) - Qq).T); keep = res < max(3.0, 3 * np.median(res[keep]))
            if keep.sum() >= 6:
                tps = RBFInterpolator(Pq[keep], Qq[keep], kernel="thin_plate_spline", smoothing=2.0)
                def warp(q):
                    ext = tps(np.array(q.exterior.coords)); ints = [tps(np.array(r.coords)) for r in q.interiors]
                    from shapely.geometry import Polygon as _P
                    w = _P(ext, ints); return w if w.is_valid else make_valid(w)
                W2 = [warp(q) for q in polys]
                cg0 = list(cur.geometry)
                def score(Wl):
                    return sum(max((g.intersection(w).area / g.union(w).area if g.intersects(w) else 0) for w in Wl) for g in cg0)
                if score(W2) > score(W_) * 1.05:
                    W_ = W2; rms = round(float(np.median(res[keep])), 2); frame_note = f"TPS on {int(keep.sum())} label pairs"
                    PAIRS[(vc, sv)] = {"sheet": Pq[keep].tolist(), "layer": Qq[keep].tolist()}
                    srow.append({"village": vc, "survey": sv, "status": frame_note})
    # one-to-one by IoU
    cg = list(cur.geometry); C = np.zeros((len(cg), len(W_))); COST = np.zeros_like(C)
    for i, g in enumerate(cg):
        for k, w in enumerate(W_):
            if Mx is None: C[i, k] = 1.0; continue
            iou = g.intersection(w).area / g.union(w).area if g.intersects(w) else 0.0
            C[i, k] = iou
            cd = g.representative_point().distance(w.representative_point())
            COST[i, k] = (1 - iou) + 0.5 * min(1.0, cd / 10.0) + 0.3 * min(2.0, abs(math.log(max(g.area, 0.1) / max(w.area, 0.1))))
    ri, ci = linear_sum_assignment(COST)
    m = {i: (k, C[i, k]) for i, k in zip(ri, ci)}
    for i, r in enumerate(cur.itertuples()):
        k, iou = m.get(i, (None, 0.0))
        printed = pn.get(k, []) if k is not None else []
        inside = sorted({norm(t) for t in printed if not t.endswith("*")}); beside = sorted({norm(t.rstrip("*")) for t in printed if t.endswith("*")})
        pr = inside or beside; lab = norm(r.lab)
        if lab and lab not in pr and lab in beside: pr = beside
        cd = cg[i].representative_point().distance(W_[k].representative_point()) if k is not None else 99
        if k is None or (iou < 0.3 and cd > 3.0): v = "UNMATCHED"
        elif not pr and not lab: v = "ok (no number)"
        elif not pr: v = "NOT-PRINTED"
        elif lab in pr: v = "ok" if len(pr) == 1 else "ok (+other number)"
        elif not lab: v = "MISSING-NUMBER"
        else: v = "WRONG"
        rows.append({"village": vc, "survey": sv, "kide": r.kide, "layer_label": r.lab, "pdf_numbers": "|".join(printed), "verdict": v, "iou": round(iou, 2),
                     "pdf_plot": k, "pdf_area": round(polys[k].area, 1) if k is not None else None, "layer_area": round(r.geometry.area, 1), "fit_rms": round(rms, 2)})
    labs_layer = {norm(x) for x in cur.lab}
    for t in sorted({norm(t) for t, p in nums} - labs_layer):
        where = [k for k, q in enumerate(polys) if q.distance([p for tt, p in nums if norm(tt) == t][0]) <= 3.0]
        rows.append({"village": vc, "survey": sv, "kide": "", "layer_label": "", "pdf_numbers": t, "verdict": "PRINTED NUMBER NOT IN LAYER", "pdf_plot": where[0] if where else None,
                     "pdf_area": round(polys[where[0]].area, 1) if where else None, "fit_rms": round(rms, 2)})
    for x in sorted(labs_layer - reg - {""}) if reg else []:
        rows.append({"village": vc, "survey": sv, "kide": f"{sv}/{x}", "layer_label": x, "verdict": "LABEL NOT IN REGISTER", "fit_rms": round(rms, 2)})
    dup = cur.lab[cur.lab != ""].value_counts()
    for x, n in dup[dup > 1].items():
        rows.append({"village": vc, "survey": sv, "kide": f"{sv}/{x}", "layer_label": x, "verdict": f"DUPLICATE LABEL x{n}", "fit_rms": round(rms, 2)})
    if os.environ.get("ANN_GPKG"):
        for k, w in enumerate(W_):
            WARPED.append({"village": vc, "survey": sv, "pdf_plot": k, "numbers": "|".join(pn[k]), "matched": k in set(ci), "geometry": w})
    used = set(ci)

    for k, q in enumerate(polys):
        prow.append({"village": vc, "survey": sv, "pdf_plot": k, "numbers": "|".join(pn[k]), "area": round(q.area, 1), "matched": k in used,
                     "x_pt": round(q.representative_point().x * PT), "y_pt": round(-q.representative_point().y * PT)})
    srow.append({"village": vc, "survey": sv, "status": "ok", "layer_plots": len(cur), "pdf_plots": len(polys), "numbers": len(nums), "fit_rms": round(rms, 2)})
    print(f"{vc} {sv}: layer {len(cur)} | PDF plots {len(polys)} | numbers {len(nums)} | fit {rms:.2f} m", flush=True)
R = pd.DataFrame(rows); R.to_csv(f"qc/{PFX}_read1.csv", index=False, encoding="utf-8-sig")
if os.environ.get("ANN_PAIRS"): json.dump({"pairs": {f"{k[0]}|{k[1]}": v for k, v in PAIRS.items()}, "rigid": {f"{k[0]}|{k[1]}": v for k, v in RIGID.items()}}, open(os.environ["ANN_PAIRS"], "w"))
if os.environ.get("ANN_GPKG"): gpd.GeoDataFrame(WARPED, geometry="geometry", crs=32644).to_file(os.environ["ANN_GPKG"], driver="GPKG")
pd.DataFrame(prow).to_csv(f"qc/{PFX}_pdfplots.csv", index=False, encoding="utf-8-sig"); pd.DataFrame(srow).to_csv(f"qc/{PFX}_sheets.csv", index=False, encoding="utf-8-sig")
print(R.verdict.value_counts().to_dict())
