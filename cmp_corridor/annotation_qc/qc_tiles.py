"""Dual label QC, step 1: for every survey of Akash's selection, draw the sheet's plot polygons (PDF line work; sketch vector when the PDF
line work does not close) on the PDF render, each with the label the LAYER carries for that plot (matched by overlap after a rigid fit of
the sheet onto the layer), as "<label> #<qc id>". Tiles -> qc/tiles/<vc>_<sv>_tNN.png; qc/plots.csv maps qc id -> layer kide.
Single-plot sheets are listed, not tiled."""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
import fiona, pandas as pd
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None
OUT = os.environ.get("TILE_OUT", "qc"); os.makedirs(f"{OUT}/tiles", exist_ok=True); os.makedirs(f"{OUT}/clean", exist_ok=True)
LAY = set(fiona.listlayers(SK))
M = gpd.read_file(os.environ.get("TILE_FILE", "D:/Projects/CMP_Roads_LPS/FMB_Georeferenced/CMP_FMB_Fine_Tunned_Merged.geojson")).to_crs(32644)
VILL = set(os.environ.get("TILE_VILLAGES", "").split(",")) - {""}
surveys = sorted({(str(v), str(k).split("/")[0]) for v, k in zip(M.village_code, M.kide) if (not VILL or str(v) in VILL) and str(k) not in ("nan", "None", "")})
M["kide"] = M.kide.astype(str); M["vc"] = M.village_code.astype(str); M["sv"] = M.kide.str.split("/").str[0]
M["lab"] = M.subdiv_no.astype(str).str.strip().replace({"None": "", "nan": "", "-": ""})
font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 15)
DPI = 170; k = DPI / 72.0
rows, single, skipped = [], [], []


def pdf_polys(pg):
    W, H = pg.rect.width, pg.rect.height
    lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2); col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1":
            continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    return [p for p in polygonize(unary_union(lines)) if p.area > 1.0]


def dens(g, step=1.5):
    g = as_poly(g); L = g.exterior; n = max(int(L.length / step), 8)
    return np.array([L.interpolate(i / n, normalized=True).coords[0] for i in range(n)])


def fit(Us, Uc, pairs):
    A = dens(Us); cA = np.array(Us.centroid.coords[0]); A0 = A - cA; B = dens(Uc); cB = np.array(Uc.centroid.coords[0]); tB = cKDTree(B)
    P0 = np.array([p[0] for p in pairs]) - cA if len(pairs) >= 2 else None; Q0 = np.array([p[1] for p in pairs]) if len(pairs) >= 2 else None

    def cost(x):
        th, tx, ty = x; ct, st = math.cos(th), math.sin(th); R = np.array([[ct, st], [-st, ct]]); P = A0 @ R + [tx, ty]
        c = (np.sqrt((tB.query(P)[0] ** 2).mean()) + np.sqrt((cKDTree(P).query(B)[0] ** 2).mean())) / 2
        if P0 is not None:
            c += np.sqrt((((P0 @ R + [tx, ty]) - Q0) ** 2).sum(1).mean())
        return c
    best = min((minimize(cost, [t, *cB], method="Nelder-Mead", options={"maxiter": 1500}) for t in np.linspace(0, 2 * math.pi, 12, endpoint=False)), key=lambda r: r.fun)
    th, tx, ty = best.x; ct, st = math.cos(th), math.sin(th)
    return [ct, -st, st, ct, tx - (ct * cA[0] - st * cA[1]), ty - (st * cA[0] + ct * cA[1])], best.fun


def one(vc, sv):
    cur = M[(M.vc == vc) & (M.sv == sv)]
    pdf = f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{vc}/{sv}.pdf"
    if not len(cur) or not os.path.exists(pdf) or f"{vc}_{sv}" not in LAY:
        skipped.append((vc, sv, "no sheet")); return
    if len(cur) == 1:
        single.append((vc, sv, cur.kide.iloc[0])); return
    pg = fitz.open(pdf)[0]
    sk = gpd.read_file(SK, layer=f"{vc}_{sv}"); sk["s"] = sk.subdiv_no.astype(str).str.strip().replace({"None": "", "nan": ""})
    raw = pdf_polys(pg)
    big, refp = max(raw, key=lambda p: p.area), max(sk.geometry, key=lambda g: g.area)
    dx, dy = refp.centroid.x - big.centroid.x, refp.centroid.y - big.centroid.y
    use_vec = abs(len(raw) - len(sk)) > max(2, 0.2 * len(sk))
    sheet = [g for g in sk.geometry] if use_vec else [translate(p, dx, dy) for p in raw]
    Us = make_valid(unary_union(sheet)); Uc = make_valid(unary_union(list(cur.geometry)))
    if Uc.area > 5 * Us.area or Us.area > 5 * Uc.area:
        skipped.append((vc, sv, f"layer {Uc.area:.0f} m2 vs sheet {Us.area:.0f} m2 - not the same parcel")); return
    pairs = []
    for s in set(sk.s) - {""}:
        a, b = sk[sk.s == s], cur[cur.lab == s]
        if len(a) == 1 and len(b) == 1:
            pairs.append((a.geometry.iloc[0].centroid.coords[0], b.geometry.iloc[0].centroid.coords[0]))
    Mx, rms = fit(Us, Uc, pairs)
    dpi = min(DPI, int(14000 / max(pg.rect.width, pg.rect.height) * 72)); k = dpi / 72.0
    pix = pg.get_pixmap(dpi=dpi); img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples); clean = img.copy(); dr = ImageDraw.Draw(img)
    q = lambda x, y: ((x - dx) * PT * k, -(y - dy) * PT * k)
    pts = []
    for n, sp in enumerate(sheet):
        wp = affine_transform(sp, Mx)
        cand = [(wp.intersection(g).area / wp.union(g).area if wp.intersects(g) else 0, kd_) for g, kd_ in zip(cur.geometry, cur.kide)]
        iou, kd_ = max(cand)
        lab = M.loc[cur.index[cur.kide == kd_][0], "lab"] if iou > 0.2 else "?"
        lab = lab or "-"
        qid = f"{sv}.{n + 1}"
        rows.append({"village": vc, "survey": sv, "qc_id": qid, "layer_kide": kd_ if iou > 0.2 else "", "layer_label": lab, "iou": round(iou, 2), "fit_rms": round(rms, 1), "frame": "vector" if use_vec else "pdf"})
        poly = sp if use_vec else translate(sp, -dx, -dy)
        poly_draw = poly if not use_vec else translate(sp, -dx, -dy)
        dr.line([q(*c) for c in translate(poly_draw, dx, dy).exterior.coords], fill=(255, 0, 0), width=2)
        c = translate(poly_draw, dx, dy).representative_point(); x, y = q(c.x, c.y)
        txt = f"{lab} #{n + 1}"
        dr.rectangle([x + 2, y - 18, x + 8 + 9 * len(txt), y + 2], fill=(255, 255, 255)); dr.text((x + 4, y - 17), txt, fill=(0, 0, 255), font=font); pts.append((x, y))
    xs = [q(*v)[0] for p in sheet for v in as_poly(p).exterior.coords]; ys = [q(*v)[1] for p in sheet for v in as_poly(p).exterior.coords]
    x0, y0, x1, y1 = int(max(min(xs) - 70, 0)), int(max(min(ys) - 70, 0)), int(min(max(xs) + 70, img.width)), int(min(max(ys) + 70, img.height))
    t = 0
    for ty in range(y0, y1, 850):
        for tx in range(x0, x1, 850):
            bx = (tx, ty, min(tx + 1000, x1), min(ty + 1000, y1))
            if [p for p in pts if bx[0] <= p[0] <= bx[2] and bx[1] <= p[1] <= bx[3]]:
                img.crop(bx).save(f"{OUT}/tiles/{vc}_{sv}_t{t:02d}.png"); clean.crop(bx).save(f"{OUT}/clean/{vc}_{sv}_t{t:02d}.png"); t += 1
    print(f"{vc} {sv}: {len(sheet)} sheet plots ({'vector' if use_vec else 'pdf'}), {len(cur)} layer plots, fit {rms:.1f} m, tiles {t}", flush=True)
for vc, sv in surveys:
    try:
        one(vc, sv)
    except Exception as e:
        skipped.append((vc, sv, "error: " + str(e)[:120])); print(vc, sv, "ERROR", e, flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/plots.csv", index=False)
json.dump({"single": single, "skipped": skipped}, open(f"{OUT}/single_skipped.json", "w"))
print("single-plot sheets:", single); print("skipped:", skipped)
print("tiles:", len(os.listdir(f"{OUT}/tiles")))
