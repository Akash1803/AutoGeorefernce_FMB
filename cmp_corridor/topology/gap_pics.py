"""Satellite picture of selected gaps (rows of a gap_report csv) with parcel outlines + labels; also which warped PDF plots/faces cover each gap.
Usage: gap_pics.py <gaps.csv> <geojson> <outprefix> <row,row,...>"""
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "D:/code/FMB_to_GeoJSON")
import geopandas as gpd, pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely import make_valid
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from PIL import Image
from autogeoref.satcheck import render_from_qgis
CSV, F, PFX, ROWS = sys.argv[1], sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4].split(",")]
D = pd.read_csv(CSV); P = gpd.read_file(F); P["kide"] = P.kide.astype(str)
G = [make_valid(g) for g in P.geometry]; U = unary_union(G)
holes = [Polygon(r) for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
W = None
for wf in ("chk_warped.gpkg", "k2_warped.gpkg"):
    if os.path.exists(wf):
        w = gpd.read_file(wf); W = w if W is None else pd.concat([W, w], ignore_index=True)
for r in ROWS:
    row = D.iloc[r]; h = min(holes, key=lambda z: abs(z.area - row.area_m2) + z.distance(Polygon([(row.x, row.y)] * 4).buffer(0.1)))
    b = h.bounds; pad = max(25, 0.6 * max(b[2] - b[0], b[3] - b[1]))
    bounds = (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)
    jpg = os.path.abspath(f"{PFX}_{r}.jpg"); render_from_qgis(bounds, jpg)
    im = Image.open(jpg); ext = [bounds[0], bounds[2], bounds[1], bounds[3]]
    fig, ax = plt.subplots(figsize=(11, 11 * im.height / im.width)); ax.imshow(im, extent=ext)
    sub = P[P.intersects(box(*bounds))]
    for g, k, st in zip(sub.geometry, sub.kide, sub.status.astype(str) if "status" in sub else [""] * len(sub)):
        col = "lime" if str(st).startswith("green") else "yellow"
        for pg in getattr(g, "geoms", [g]):
            x, y = pg.exterior.xy; ax.plot(x, y, color=col, lw=1.1)
        c = g.representative_point(); ax.text(c.x, c.y, k.split("/")[-1] if "/" in k else k, color="white", fontsize=7, ha="center", va="center", bbox=dict(fc="black", alpha=0.35, lw=0, pad=0.5))
    x, y = h.exterior.xy; ax.fill(x, y, color="red", alpha=0.35); ax.plot(x, y, color="red", lw=1.5)
    if os.environ.get("GP_ADDS"):
        from shapely import wkt as swkt
        for a in json.load(open(os.environ["GP_ADDS"]))["add"]:
            ga = swkt.loads(a["wkt"])
            if ga.intersects(box(*bounds)):
                for pg in getattr(ga, "geoms", [ga]): xx, yy = pg.exterior.xy; ax.fill(xx, yy, color="cyan", alpha=0.45); ax.plot(xx, yy, color="cyan", lw=1.2)
    cover = ""
    if W is not None:
        wi = W[W.intersects(h)]
        parts = [f"{v[-3:]}:{s} plot {p} [{n}] {100 * g.intersection(h).area / h.area:.0f}%" for v, s, p, n, g in zip(wi.village, wi.survey, wi.pdf_plot, wi.numbers.astype(str), wi.geometry) if g.intersection(h).area > 0.5]
        cover = " | ".join(parts)
    ax.set_title(f"gap {r}: {row.area_m2} m2, width {row.width_m} m, {row.in_buffer_pct:.0f}% in buffer\nPDF cover: {cover[:160]}", fontsize=9)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_axis_off(); fig.tight_layout(); fig.savefig(f"{PFX}_{r}.png", dpi=110); plt.close(fig)
    print(r, round(h.area, 1), "cover:", cover)
