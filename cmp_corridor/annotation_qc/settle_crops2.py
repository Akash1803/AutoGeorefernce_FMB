"""Settlement crops v2: per non-OK, non-green, not-yet-decided decision row, the PDF around the MATCHED PDF plot (drawn red from the
sheet's own line work, no warp) or, without a match, the layer plot found by (village, kide, area) inverse-mapped. 6 per image."""
import sys, os, json, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd, geopandas as gpd, fitz
from scipy.interpolate import RBFInterpolator
from shapely.geometry import LineString
from shapely.ops import unary_union, polygonize
from PIL import Image, ImageDraw, ImageFont
PT = 72 / 25.4 * 2
D = pd.read_csv("qc/two_decisions.csv", dtype=str).fillna("")
R1 = pd.read_csv("qc/two_read1.csv", dtype=str).fillna("")
done = set(pd.read_csv("qc/settle/decided_keys.csv", dtype=str).key)
A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); A["key"] = A.village_code.astype(str) + "|" + A.kide.astype(str)
green = set(A.key[A.status.astype(str).str.startswith("green")])
D = D[~D.status.isin(["OK"]) & ~D.key.isin(green) & ~D.key.isin(done)].copy()
PAIRS = json.load(open("chk_pairs.json")); L = A.to_crs(32644) if A.crs and A.crs.to_epsg() != 32644 else A
os.makedirs("qc/settle2", exist_ok=True)
font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 16); small = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 13)
cache = {}
def polys_of(vc, sv):
    if (vc, sv) in cache: return cache[(vc, sv)]
    pg = fitz.open(f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{vc}/{sv}.pdf")[0]; W, H = pg.rect.width, pg.rect.height; lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2); col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1": continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    cache[(vc, sv)] = (pg, [p for p in polygonize(unary_union(lines)) if p.area > 0.5]); return cache[(vc, sv)]
def inverse_map(vc, sv):
    k = f"{vc}|{sv}"
    if k in PAIRS["pairs"]:
        p = PAIRS["pairs"][k]; tps = RBFInterpolator(np.array(p["layer"]), np.array(p["sheet"]), kernel="thin_plate_spline", smoothing=2.0); return lambda xy: tps(np.asarray(xy))
    if k in PAIRS["rigid"]:
        a, b, c, d, e, f = PAIRS["rigid"][k]; Ai = np.linalg.inv(np.array([[a, b], [c, d]])); t = np.array([e, f]); return lambda xy: (np.asarray(xy) - t) @ Ai.T
    return None
panels, index = [], []
for r in D.itertuples():
    if not os.path.exists(f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{r.village}/{r.survey}.pdf"): index.append({"key": r.key, "image": "", "panel": "", "why": "no PDF"}); continue
    pg, polys = polys_of(r.village, r.survey); ring = None
    if r.pdf_plot not in ("", "nan"):
        k = int(float(r.pdf_plot))
        if k < len(polys): ring = np.array(polys[k].exterior.coords)
    if ring is None:
        inv = inverse_map(r.village, r.survey)
        rr = R1[(R1.village == r.village) & (R1.kide == r.kide)]
        g = L[L.key == r.key]
        if inv is not None and len(g):
            # pick the layer row whose area matches the read1 row
            tgt = float(rr.layer_area.iloc[0]) if len(rr) and rr.layer_area.iloc[0] else None
            geo = min(g.geometry, key=lambda x: abs(x.area - tgt)) if tgt else max(g.geometry, key=lambda x: x.area)
            ring = inv(np.array(geo.exterior.coords))
    if ring is None: index.append({"key": r.key, "image": "", "panel": "", "why": "no pose"}); continue
    cx, cy = ring[:, 0].mean() * PT, -ring[:, 1].mean() * PT; half = max(60, 0.7 * max(np.ptp(ring[:, 0]), np.ptp(ring[:, 1])) * PT)
    panels.append((r, pg, ring, cx, cy, half))
Z = 3
for i in range(0, len(panels), 6):
    group = panels[i:i + 6]; tiles = []
    for r, pg, ring, cx, cy, half in group:
        clip = fitz.Rect(cx - half, cy - half, cx + half, cy + half) & pg.rect
        pix = pg.get_pixmap(dpi=72 * Z, clip=clip); im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).resize((480, 480))
        sx, sy = 480 / clip.width, 480 / clip.height; dr = ImageDraw.Draw(im)
        pts = [((x * PT - clip.x0) * sx, (-y * PT - clip.y0) * sy) for x, y in ring]; dr.line(pts + [pts[0]], fill=(255, 0, 0), width=3)
        dr.rectangle([0, 0, 480, 40], fill=(255, 255, 255)); dr.text((4, 2), f"{r.village[-3:]} {r.kide}  layer={r.layer_label or '(none)'}", fill=(200, 0, 0), font=font); dr.text((4, 22), f"M={r.machine} A={r.A} B={r.B}  [{r.status[:12]}]", fill=(0, 0, 160), font=small)
        tiles.append(im)
    Wd = Image.new("RGB", (3 * 484, 2 * 484), (40, 40, 40))
    for j, im in enumerate(tiles): Wd.paste(im, ((j % 3) * 484, (j // 3) * 484))
    name = f"qc/settle2/s_{i // 6:02d}.png"; Wd.save(name)
    for j, (r, *_x) in enumerate(group): index.append({"key": r.key, "image": name, "panel": j + 1, "status": r.status, "layer_label": r.layer_label, "machine": r.machine, "A": r.A, "B": r.B, "new_label": r.new_label})
pd.DataFrame(index).to_csv("qc/settle2/index.csv", index=False, encoding="utf-8-sig")
print("panels:", len(panels), "images:", (len(panels) + 5) // 6, "| no pose:", sum(1 for x in index if x.get("image") == ""))
