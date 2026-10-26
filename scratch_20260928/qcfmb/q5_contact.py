"""Dual QC, step 5: clean crops of the PDF for the rows that need a look, 6 per image; the plot is outlined thinly in red,
nothing is printed over the sheet, the caption carries only the row id. -> out/cs/v_<k>.png, out/cs_index.csv"""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, pymupdf, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import affinity
from shapely.ops import unary_union
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
V = pd.read_csv(os.path.join(HERE, "out", sys.argv[1] if len(sys.argv) > 1 else "visual.csv"), dtype=str).fillna(""); tag = sys.argv[2] if len(sys.argv) > 2 else "v"
L = pd.read_csv(os.path.join(HERE, "out", "labels5_all.csv"), dtype=str).fillna("")
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
cache = {}
def sheet(v, s):
    if (v, s) not in cache:
        if len(cache) > 12: cache.clear()
        M = scales.get((v, s), 500) * 0.0254 / 72.0; page = pymupdf.open(os.path.join(R, "FMB_Sketches", v, s + ".pdf"))[0]; H = page.rect.height
        S = gpd.read_file(os.path.join(R, "FMB_Vector", v, s + "_parcels.geojson")); S["geometry"] = S.geometry.buffer(0)
        cache[(v, s)] = (page, H, {str(int(p)): affinity.affine_transform(g, [1 / M, 0, 0, -1 / M, 0, H]) for p, g in zip(S.poly_id, S.geometry)})
    return cache[(v, s)]
V = V.sort_values(["village_code", "sheet", "row"]).reset_index(drop=True); os.makedirs(os.path.join(HERE, "out", "cs"), exist_ok=True); idx = []
for k in range(0, len(V), 6):
    chunk = V.iloc[k:k + 6]; fig, axs = plt.subplots(2, 3, figsize=(16.5, 11.6)); axs = np.array(axs).ravel()
    for ax in axs: ax.axis("off")
    for ax, r in zip(axs, chunk.itertuples()):
        page, H, P = sheet(r.village_code, r.sheet); pids = [p for p in r.pids.split("/") if p in P]
        if not pids: ax.set_title(f"#{r.row} (no polygon)", fontsize=13); continue
        g = unary_union([P[p] for p in pids]); b = g.bounds
        c = L[(L.village_code == r.village_code) & (L.sheet == r.sheet) & (L.plotno == "True") & (L.pid.isin(pids) | L.near_pid.isin(pids))]
        hw = max(b[2] - b[0], b[3] - b[1]) / 2 + 14
        if hw <= 120 or not len(c): cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2; hw = min(hw, 160)
        else: cx, cy = float(c.px.astype(float).mean()), float(c.py.astype(float).mean()); hw = 110                 # a big plot: look where the numbers stand
        hw = max(hw, 42); clip = pymupdf.Rect(cx - hw, cy - hw, cx + hw, cy + hw) & page.rect
        pix = page.get_pixmap(dpi=int(min(560 / (2 * hw) * 72, 600)), clip=clip); img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
        ax.imshow(img, extent=(clip.x0, clip.x1, clip.y1, clip.y0))
        for q in (g.geoms if g.geom_type == "MultiPolygon" else [g]):
            xs, ys = np.array(q.exterior.coords).T; ax.plot(xs, ys, c="red", lw=2.4, alpha=0.3)
        ax.set_xlim(clip.x0, clip.x1); ax.set_ylim(clip.y1, clip.y0); ax.set_title(f"#{r.row}", fontsize=14, color="darkred")
        idx.append({"image": f"{tag}_{k // 6:03d}.png", "row": r.row, "village_code": r.village_code, "sheet": r.sheet, "pids": r.pids})
    plt.tight_layout(); plt.savefig(os.path.join(HERE, "out", "cs", f"{tag}_{k // 6:03d}.png"), dpi=62); plt.close()
pd.DataFrame(idx).to_csv(os.path.join(HERE, "out", f"cs_index_{tag}.csv"), index=False); print("rows", len(V), "| images", (len(V) + 5) // 6)
