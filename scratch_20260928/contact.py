"""Contact sheets: one clean crop of the PDF per plot (outline only, nothing drawn over the print), caption above.
usage: python contact.py <name> <rid,rid,...> [cols]   -> out/cs/<name>_<k>.png"""
import warnings, sys, os, json; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, pymupdf, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import affinity
from autogeoref import visible, fabric_check as fc
name, rids = sys.argv[1], [r for r in sys.argv[2].split(",") if r]; cols = int(sys.argv[3]) if len(sys.argv) > 3 else 3
R = r"D:\Projects\Tambaram_Chengalpattu"; D = pd.read_csv("out/decide.csv", dtype=str).fillna("").set_index("rid"); AL = json.load(open("out/align2.json"))
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
cache = {}
def sheet(V, s):
    if (V, s) not in cache:
        M = scales.get((V, s), 500) * 0.0254 / 72.0; page = pymupdf.open(os.path.join(R, "FMB_Sketches", V, s + ".pdf"))[0]; H = page.rect.height
        S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); S["geometry"] = S.geometry.buffer(0); sh = visible.shown(V)
        L = gpd.read_file(sh[s][0], layer=sh[s][1]); L = L.to_crs(32644) if L.crs.to_epsg() != 32644 else L; L = L[L.geometry.notna() & ~L.geometry.is_empty]
        cache[(V, s)] = (M, page, H, {int(p): g for p, g in zip(S.poly_id, S.geometry)}, L)
    return cache[(V, s)]
per = cols * 2; os.makedirs("out/cs", exist_ok=True)
for k in range(0, len(rids), per):
    chunk = rids[k:k + per]; fig, axs = plt.subplots(2, cols, figsize=(5.2 * cols, 11)); axs = np.array(axs).ravel()
    for ax in axs: ax.axis("off")
    for ax, rid in zip(axs, chunk):
        r = D.loc[rid]; V, s = r.village_code, r.survey; M, page, H, SP, L = sheet(V, s); topdf = lambda g: affinity.affine_transform(g, [1 / M, 0, 0, -1 / M, 0, H])
        if r.pid_ok == "True" and r.pid and int(r.pid) in SP: g = topdf(SP[int(r.pid)]); how = ""
        else:
            a = AL[f"{V}|{s}"]; Rm = np.array([[a[0], a[1]], [a[2], a[3]]]); t = np.array([a[4], a[5]]); Ri = Rm.T; ti = -Ri @ t
            lg = max(fc._polys(L.geometry.loc[int(r.fid)].buffer(0)), key=lambda q: q.area); g = topdf(affinity.affine_transform(lg, [Ri[0, 0], Ri[0, 1], Ri[1, 0], Ri[1, 1], ti[0], ti[1]])); how = " ~placed"
        b = g.bounds; w = (max(b[2] - b[0], b[3] - b[1], 60) / 2 + 22) * float(os.environ.get('WIDE', '1')); cx, cy = (b[0] + b[2]) / 2, (b[1] + b[3]) / 2
        pix = page.get_pixmap(dpi=int(min(520 / (2 * w) * 72, 700)), clip=pymupdf.Rect(cx - w, cy - w, cx + w, cy + w)); img = np.frombuffer(pix.samples, np.uint8).reshape(pix.h, pix.w, pix.n)[:, :, :3]
        ax.imshow(img, extent=(cx - w, cx + w, cy + w, cy - w)); xs, ys = np.array(g.exterior.coords).T; ax.plot(xs, ys, c="red", lw=2.2, alpha=0.35)
        ax.set_title(f"#{rid} {s}: layer [{r.layer_plot or ' '}] -> {r.proposed or '-'}  ({float(r.area):.0f} m2){how}", fontsize=11, color="darkred")
    out = f"out/cs/{name}_{k // per}.png"; plt.tight_layout(); plt.savefig(out, dpi=72); plt.close(); print(out, chunk)
