"""Clean view of a sheet: the PDF page (or a crop around sheet polygons), the chosen polygons outlined thinly.
usage: python sheet_view.py <village_code> <survey> <pid,pid,...|all> <out.png> [margin_pt] [dpi]"""
import warnings, sys, os; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, pymupdf, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import affinity
from shapely.ops import unary_union
V, s, pids, out = sys.argv[1:5]; mar = float(sys.argv[5]) if len(sys.argv) > 5 else 30; dpi = int(sys.argv[6]) if len(sys.argv) > 6 else 200
R = r"D:\Projects\Tambaram_Chengalpattu"
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
M = scales.get((V, s), 500) * 0.0254 / 72.0; page = pymupdf.open(os.path.join(R, "FMB_Sketches", V, s + ".pdf"))[0]; H = page.rect.height
S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); S["geometry"] = S.geometry.buffer(0)
topdf = lambda g: affinity.affine_transform(g, [1 / M, 0, 0, -1 / M, 0, H])
sel = S if pids == "all" else S[S.poly_id.astype(int).isin([int(p) for p in pids.split(",")])]
G = [topdf(g) for g in sel.geometry]; b = unary_union(G).bounds
clip = pymupdf.Rect(max(b[0] - mar, 0), max(b[1] - mar, 0), min(b[2] + mar, page.rect.width), min(b[3] + mar, H))
pix = page.get_pixmap(dpi=dpi, clip=clip); img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
fig = plt.figure(figsize=(pix.w / 100, pix.h / 100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.imshow(img, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
if pids != "all":
    for g, p in zip(G, sel.poly_id):
        xs, ys = np.array(g.exterior.coords).T; ax.plot(xs, ys, c="red", lw=2.5, alpha=0.3)
ax.set_xlim(clip.x0, clip.x1); ax.set_ylim(clip.y1, clip.y0); plt.savefig(out, dpi=100); print(out, pix.w, pix.h, [(int(p), round(a, 1)) for p, a in zip(sel.poly_id, sel.geometry.area)])
