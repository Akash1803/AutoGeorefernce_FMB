"""Draw every placed plot (as QGIS shows it) back onto its own PDF sheet with its current plot number, so the printed
number can be read beside it. Red text = current layer value ('?' = blank), '#k' = row id in out/xcheck.csv.
usage: python overlay_check.py <village_code> <survey> [tile]    -> out/ov/<vvv>_<survey>[_tile].png"""
import warnings, sys, os; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely, pymupdf, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import affinity
from shapely.ops import unary_union
from autogeoref import visible, fabric_check as fc
V, s = sys.argv[1], sys.argv[2]; tile = sys.argv[3] if len(sys.argv) > 3 else ""
R = r"D:\Projects\Tambaram_Chengalpattu"
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
sc = scales.get((V, s), 500); M = sc * 0.0254 / 72.0
page = pymupdf.open(os.path.join(R, "FMB_Sketches", V, s + ".pdf"))[0]; H = page.rect.height
S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); S["geometry"] = S.geometry.buffer(0)
sh = visible.shown(V); L = gpd.read_file(sh[s][0], layer=sh[s][1]); L = L.to_crs(32644) if L.crs.to_epsg() != 32644 else L; L = L[L.geometry.notna() & ~L.geometry.is_empty]
X = pd.read_csv("out/xcheck2.csv", dtype=str).fillna(""); X = X[(X.village_code == V) & (X.survey == s)]; rid = {int(f): i for i, f in zip(X.index, X.fid)}
ocr = {int(f): o for f, o in zip(X.fid, X.sheet_plot)}; inb = {int(f): b == "True" for f, b in zip(X.fid, X.in_buffer)}
ONLY = set(int(k) for k in os.environ.get("ONLY", "").split(",") if k)
body_s = unary_union(list(S.geometry)); body_l = unary_union([g.buffer(0) for g in L.geometry])
_, (Rm, t) = fc.shape_iou(body_s, body_l)
import json as _j
_a = _j.load(open("out/align2.json")).get(f"{V}|{s}")
if _a and os.environ.get("ALIGN2", "1") == "1": Rm = np.array([[_a[0], _a[1]], [_a[2], _a[3]]]); t = np.array([_a[4], _a[5]])
Ri = Rm.T; ti = -Ri @ t
inv = [Ri[0, 0], Ri[0, 1], Ri[1, 0], Ri[1, 1], ti[0], ti[1]]
topdf = lambda g: affinity.affine_transform(affinity.affine_transform(g, inv), [1 / M, 0, 0, -1 / M, 0, H])
P = [(fid, topdf(r.geometry.buffer(0)), r.get("plot_no")) for fid, r in L.iterrows()]
if os.environ.get("BY_PID"):                                   # draw the sheet's own polygons, labelled with the layer plot of the same poly_id
    toS = lambda g: affinity.affine_transform(g, [1 / M, 0, 0, -1 / M, 0, H]); sp = {int(q.poly_id): q.geometry for _, q in S.iterrows()}
    P = [(fid, toS(sp[int(r.poly_id)]), r.get("plot_no")) for fid, r in L.iterrows() if int(r.poly_id) in sp]
b = np.array(unary_union([p for _, p, _ in P] + [affinity.affine_transform(body_s, [1 / M, 0, 0, -1 / M, 0, H])]).bounds)
x0, y0, x1, y1 = b[0] - 25, b[1] - 25, b[2] + 25, b[3] + 25
if tile.startswith("z"):                                   # z<rid>[,<rid>...] : zoom on these rows
    ids = [int(k) for k in tile[1:].split(",")] if tile != "zbuf" else [rid[f] for f in rid if inb.get(f)]; zb = np.array(unary_union([p for fid, p, _ in P if rid.get(fid) in ids]).bounds); PAD = float(os.environ.get('PAD', '40')); x0, y0, x1, y1 = zb[0] - PAD, zb[1] - PAD, zb[2] + PAD, zb[3] + PAD
elif tile:
    nx, ny, ix, iy = map(int, tile.split(",")); w, h = (x1 - x0) / nx, (y1 - y0) / ny; x0, y0 = x0 + ix * w, y0 + iy * h; x1, y1 = x0 + w, y0 + h
dpi = min(2600 / max(x1 - x0, y1 - y0) * 72, 600)
pix = page.get_pixmap(dpi=int(dpi), clip=pymupdf.Rect(x0, y0, x1, y1)); img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
fig, ax = plt.subplots(figsize=(20, 20 * (y1 - y0) / (x1 - x0)))
ax.imshow(img, extent=(x0, x1, y1, y0))
for fid, p, pn in ([] if os.environ.get('NODRAW') else P):
    for q in getattr(p, "geoms", [p]):
        xs, ys = np.array(q.exterior.coords).T; ax.plot(xs, ys, c="red", lw=0.5, alpha=0.5)
    c = p.representative_point(); lab = "?" if pn is None or str(pn) in ("nan", "None", "") else str(pn)
    if os.environ.get("ONLY_BUF") and not inb.get(fid): continue
    if ONLY and rid.get(fid) not in ONLY: continue
    if x0 < c.x < x1 and y0 < c.y < y1:
        o = ocr.get(fid, ""); txt = f"{lab} #{rid.get(fid, '')}" + ("" if o == ("" if lab == "?" else lab) else f" [ocr {o or '-'}]") + (" B" if inb.get(fid) else "")
        ax.annotate(txt, (c.x, c.y), xytext=(0, -20), textcoords="offset points", color="red", fontsize=int(os.environ.get("FS", "15")), weight="bold", ha="center",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.55))
ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis("off")
os.makedirs("out/ov", exist_ok=True); out = f"out/ov/{'b_' if os.environ.get('ONLY_BUF') else ''}{V[-3:]}_{s}{'_' + tile.replace(',', '') if tile else ''}.png"
plt.savefig(out, dpi=72, bbox_inches="tight"); print(out, "| plots", len(P))
