"""qc tile ids -> PDF plot centroids (page points), reproducing qc_tiles.py's polygon list and order per survey, so that tile reads
(by qc_id) can be joined to the machine read's PDF plots (ann_*_pdfplots.csv: pdf_plot, x_pt, y_pt) by nearest centroid.
Usage: python tile_poly_map.py <tiledir>  -> <tiledir>/qc_centroids.csv"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace"); sys.path.insert(0, ".")
import pandas as pd, geopandas as gpd, fitz, fiona
from shapely.geometry import LineString
from shapely.ops import unary_union, polygonize
PT = 72 / 25.4 * 2
SK = "D:/Projects/CMP_Roads_LPS/FMB_Vector/CMP_FMB_Vector_all_sketches.gpkg"
td = sys.argv[1]; Q = pd.read_csv(f"{td}/plots.csv", dtype=str).fillna("")
LAY = set(fiona.listlayers(SK)); rows = []
for (vc, sv), grp in Q.groupby(["village", "survey"]):
    pdf = f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{vc}/{sv}.pdf"
    if not os.path.exists(pdf) or f"{vc}_{sv}" not in LAY: continue
    pg = fitz.open(pdf)[0]; W, H = pg.rect.width, pg.rect.height; lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2); col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1": continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    raw = [p for p in polygonize(unary_union(lines)) if p.area > 1.0]
    frame = grp.frame.iloc[0]
    if frame == "vector":
        sk = gpd.read_file(SK, layer=f"{vc}_{sv}")
        big, refp = max(raw, key=lambda p: p.area), max(sk.geometry, key=lambda g: g.area)
        dx, dy = refp.centroid.x - big.centroid.x, refp.centroid.y - big.centroid.y
        sheet = [g for g in sk.geometry]          # sketch frame: back to the PDF frame with -dx, -dy
        for n, sp in enumerate(sheet):
            c = sp.representative_point(); rows.append({"village": vc, "survey": sv, "qc_id": f"{sv}.{n + 1}", "x_pt": round((c.x - dx) * PT), "y_pt": round(-(c.y - dy) * PT)})
    else:
        for n, sp in enumerate(raw):
            c = sp.representative_point(); rows.append({"village": vc, "survey": sv, "qc_id": f"{sv}.{n + 1}", "x_pt": round(c.x * PT), "y_pt": round(-c.y * PT)})
pd.DataFrame(rows).to_csv(f"{td}/qc_centroids.csv", index=False); print("qc ids mapped:", len(rows), "of", len(Q))
