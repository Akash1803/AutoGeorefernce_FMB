"""Upright label-size characters inside one sheet polygon, with the distance to every candidate character."""
import sys, os, math, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read3 as R3, fmb_to_geojson as F, pdfplumber, pandas as pd, geopandas as gpd, shapely
D = pd.read_csv("out/decide.csv", dtype=str).fillna("").set_index("rid")
scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R3.R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
for rid in sys.argv[1].split(","):
    r = D.loc[rid]; V, s = r.village_code, r.survey; M = F.PT_MM * scales.get((V, s), 500) / 1000.0
    S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); g = S[S.poly_id == int(r.pid)].geometry.iloc[0].buffer(0)
    page = pdfplumber.open(os.path.join(R3.R, "FMB_Sketches", V, s + ".pdf")).pages[0]; H = float(page.height)
    hL = R3.cap_height(page, set(D[(D.village_code == V) & (D.survey == s)].layer_plot) - {""}); C = R3.curves(page); out = []
    for i, o in enumerate(C):
        cx, cy = (o["x0"] + o["x1"]) / 2, (o["y0"] + o["y1"]) / 2
        if not g.contains(shapely.Point(cx * M, (H - cy) * M)): continue
        h = o["y1"] - o["y0"]
        if abs(h - hL) > 0.06 * hL: continue
        poly = shapely.Polygon(o["pts"]).buffer(0)
        parts = [o] + [q for j, q in enumerate(C) if j != i and (q["y1"] - q["y0"]) <= 0.93 * hL and q["x0"] > o["x0"] and q["x1"] < o["x1"] and q["y0"] > o["y0"] and q["y1"] < o["y1"]
                       and shapely.Polygon(q["pts"]).buffer(0).intersection(poly).area > 0.97 * max(shapely.Polygon(q["pts"]).buffer(0).area, 1e-9)]
        bm, ar = F._bitmap(parts); cand = {}
        for ch, b_, a_ in R3.LIB + R3.NARROW:
            if abs(math.log(ar / a_)) > 0.30: continue
            d = R3.ham(bm, b_); cand[ch] = min(cand.get(ch, 999), d)
        out.append((round(cx, 1), round(cy, 1), round(h / hL, 3), round(ar, 3), len(parts) - 1, sorted(cand.items(), key=lambda t: t[1])[:4]))
    print(f"#{rid} {s} layer [{r.layer_plot}] cap {hL:.2f} pt"); [print("   ", o) for o in sorted(out)]
