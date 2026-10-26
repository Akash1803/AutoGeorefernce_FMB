"""Whole sheet or a crop, chosen sheet polygons tinted and tagged p<id>/<area> (position check, not for reading labels).
usage: python sheet_mark.py <village_code> <survey> <pid,pid,..> <out.png> [dpi] [crop_margin_pt|page]"""
import warnings, os, sys; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, pymupdf, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import affinity
from shapely.ops import unary_union
V, s, pids, out = sys.argv[1:5]; dpi = int(sys.argv[5]) if len(sys.argv) > 5 else 110; crop = sys.argv[6] if len(sys.argv) > 6 else "page"; pids = [int(p) for p in pids.split(",")]
R = r"D:\Projects\Tambaram_Chengalpattu"
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
M = scales.get((V, s), 500) * 0.0254 / 72.0; page = pymupdf.open(os.path.join(R, "FMB_Sketches", V, s + ".pdf"))[0]; H = page.rect.height
S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); S["geometry"] = S.geometry.buffer(0); S["q"] = [affinity.affine_transform(g, [1 / M, 0, 0, -1 / M, 0, H]) for g in S.geometry]
print(S[["poly_id", "plot_no"]].assign(a=S.area.round(1)).to_string(), "\nheader:", [l for l in page.get_text().splitlines() if "Area" in l or "Survey No" in l])
if crop == "page": b = unary_union(list(S.q)).bounds; mar = 40
else: b = unary_union([q for p, q in zip(S.poly_id, S.q) if int(p) in pids]).bounds; mar = float(crop)
clip = pymupdf.Rect(max(b[0] - mar, 0), max(b[1] - mar, 0), min(b[2] + mar, page.rect.width), min(b[3] + mar, H))
pix = page.get_pixmap(dpi=dpi, clip=clip); img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
fig = plt.figure(figsize=(pix.w / 100, pix.h / 100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.imshow(img, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
for p, q, g in zip(S.poly_id, S.q, S.geometry):
    if int(p) in pids:
        xs, ys = np.array(q.exterior.coords).T; ax.fill(xs, ys, fc="yellow", alpha=0.25, ec="red", lw=1.2); c = q.representative_point(); ax.text(c.x, c.y, f"p{p}\n{g.area:.0f}", color="darkgreen", fontsize=9, ha="center", weight="bold")
ax.set_xlim(clip.x0, clip.x1); ax.set_ylim(clip.y1, clip.y0); plt.savefig(out, dpi=100); print(out, pix.w, pix.h, "scale", scales.get((V, s), 500))
