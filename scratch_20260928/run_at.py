"""All blue characters on one text line of a sheet: usage run_at.py <village> <survey> <y_pt> <x0_pt> <x1_pt>"""
import sys, os, math, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read3 as R3, fmb_to_geojson as F, pdfplumber, shapely
V, s, y, x0, x1 = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
page = pdfplumber.open(os.path.join(R3.R, "FMB_Sketches", V, s + ".pdf")).pages[0]; C = R3.curves(page); out = []
ALL = R3.LIB + R3.NARROW
for i, o in enumerate(C):
    cx = (o["x0"] + o["x1"]) / 2; h = o["y1"] - o["y0"]
    if not (x0 <= cx <= x1) or h < 5 or abs(o["y0"] - (y - 3.65)) > 1.2: continue          # same cap line (top of the letters)
    poly = shapely.Polygon(o["pts"]).buffer(0)
    parts = [o] + [q for j, q in enumerate(C) if j != i and q["x0"] > o["x0"] and q["x1"] < o["x1"] and q["y0"] > o["y0"] and q["y1"] < o["y1"]
                   and shapely.Polygon(q["pts"]).buffer(0).intersection(poly).area > 0.97 * max(shapely.Polygon(q["pts"]).buffer(0).area, 1e-9)]
    bm, ar = F._bitmap(parts); cand = {}
    for ch, b_, a_ in ALL:
        if abs(math.log(ar / a_)) > 0.35: continue
        cand[ch] = min(cand.get(ch, 999), R3.ham(bm, b_))
    rk = sorted(cand.items(), key=lambda t: t[1])[:3]; out.append((round(cx, 1), round(h, 2), round(ar, 2), rk))
out.sort(); print("".join(o[3][0][0] if o[3] else "?" for o in out)); [print("  ", o) for o in out]
