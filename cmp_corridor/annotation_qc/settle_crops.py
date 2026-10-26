"""Settlement crops for every non-OK decision: the PDF around the plot, with the LAYER plot's outline drawn in red (layer -> sheet
by the inverse of the survey's label-anchored warp, or the inverse rigid pose), the layer label and the three reads in the title.
6 crops per sheet image -> qc/settle/settle_NN.png + qc/settle/index.csv (one row per crop: key, image, panel)."""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd, geopandas as gpd, fitz
from scipy.interpolate import RBFInterpolator
from shapely.affinity import affine_transform
from PIL import Image, ImageDraw, ImageFont
PT = 72 / 25.4 * 2
D = pd.read_csv("qc/two_decisions.csv", dtype=str).fillna("")
D = D[~D.status.isin(["OK"])].copy()
PAIRS = json.load(open("chk_pairs.json"))
L = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); L["kide"] = L.kide.astype(str); L["vc"] = L.village_code.astype(str)
os.makedirs("qc/settle", exist_ok=True)
font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 16); small = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 13)
def inverse_map(vc, sv):
    k = f"{vc}|{sv}"
    if k in PAIRS["pairs"]:
        p = PAIRS["pairs"][k]; tps = RBFInterpolator(np.array(p["layer"]), np.array(p["sheet"]), kernel="thin_plate_spline", smoothing=2.0)
        return lambda xy: tps(np.asarray(xy))
    if k in PAIRS["rigid"]:
        a, b, c, d, e, f = PAIRS["rigid"][k]; A = np.array([[a, b], [c, d]]); t = np.array([e, f]); Ai = np.linalg.inv(A)
        return lambda xy: (np.asarray(xy) - t) @ Ai.T
    return None
panels, index = [], []
for r in D.itertuples():
    inv = inverse_map(r.village, r.survey)
    g = L[(L.vc == r.village) & (L.kide == r.kide)]
    if inv is None or not len(g): index.append({"key": r.key, "image": "", "panel": "", "why": "no pose"}); continue
    geo = g.geometry.iloc[0] if len(g) == 1 else max(g.geometry, key=lambda x: x.area)
    ring = inv(np.array(geo.exterior.coords))
    cx, cy = ring[:, 0].mean() * PT, -ring[:, 1].mean() * PT
    half = max(60, 0.7 * max(np.ptp(ring[:, 0]), np.ptp(ring[:, 1])) * PT)
    panels.append((r, ring, cx, cy, half))
Z = 3
for i in range(0, len(panels), 6):
    group = panels[i:i + 6]; tiles = []
    for r, ring, cx, cy, half in group:
        pg = fitz.open(f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{r.village}/{r.survey}.pdf")[0]
        clip = fitz.Rect(cx - half, cy - half, cx + half, cy + half) & pg.rect
        pix = pg.get_pixmap(dpi=72 * Z, clip=clip); im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).resize((480, 480))
        sx, sy = 480 / clip.width, 480 / clip.height; dr = ImageDraw.Draw(im)
        pts = [((x * PT - clip.x0) * sx, (-y * PT - clip.y0) * sy) for x, y in ring]
        dr.line(pts + [pts[0]], fill=(255, 0, 0), width=3)
        title = f"{r.village[-3:]} {r.kide}  layer={r.layer_label or '(none)'}"; sub = f"M={r.machine} A={r.A} B={r.B}"
        dr.rectangle([0, 0, 480, 40], fill=(255, 255, 255)); dr.text((4, 2), title, fill=(200, 0, 0), font=font); dr.text((4, 22), sub, fill=(0, 0, 160), font=small)
        tiles.append(im)
    W = Image.new("RGB", (3 * 484, 2 * 484), (40, 40, 40))
    for j, im in enumerate(tiles): W.paste(im, ((j % 3) * 484, (j // 3) * 484))
    name = f"qc/settle/settle_{i // 6:02d}.png"; W.save(name)
    for j, (r, *_rest) in enumerate(group): index.append({"key": r.key, "image": name, "panel": j + 1, "status": r.status, "layer_label": r.layer_label, "machine": r.machine, "A": r.A, "B": r.B})
pd.DataFrame(index).to_csv("qc/settle/index.csv", index=False, encoding="utf-8-sig")
print("crops:", len(panels), "images:", (len(panels) + 5) // 6, "| without pose:", sum(1 for x in index if x.get("image") == ""))
