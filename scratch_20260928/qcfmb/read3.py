"""Third read of the printed plot numbers, built for labels that OVERLAP other text.
On these sheets every character is its own filled outline, so overlapping text never changes a character's shape - it
only confuses which characters belong together. So: classify each character by its own outline (+ its true inner
counters), keep upright characters of the label font size only (bounding-box height = the cap height of the labels the
first read already confirmed), and join characters that share a baseline. Alphabet = digits + A-Z (sheet font outlines).
-> out/labels3.csv"""
import sys, os, json, math, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely, pdfplumber
import fmb_to_geojson as F
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
FONT = json.load(open(os.path.join(HERE, "out", "font_bitmaps.json")))
BASE = [(c, tuple(b), a) for c, b, a in json.load(open(r"D:\code\FMB_to_GeoJSON\glyphlib.json")) if c not in "~()"]
LIB = BASE + [(c, tuple(FONT[c][0]), FONT[c][1]) for c in "0123456789ABCDEFGHKLMNPRSTUVWXYZ"]
NARROW = [(c, tuple(FONT[c][0]), FONT[c][1]) for c in "IJ"]
ham = lambda a, b: sum(1 for x, y in zip(a, b) if x != y)


def curves(page):
    out = []
    for c in page.curves:
        if c.get("non_stroking_color") != (0.0, 0.0, 1.0): continue
        pts = F.flatten(c["path"])
        if len(pts) < 3: continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        out.append({"pts": pts, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys)})
    for r in page.rects:                      # the letter I is a plain rectangle: pdfplumber files it under rects, not curves
        if r.get("non_stroking_color") != (0.0, 0.0, 1.0) or not r.get("fill", True): continue
        x0, x1, y0, y1 = r["x0"], r["x1"], r["top"], r["bottom"]
        if x1 - x0 < 0.3 or y1 - y0 < 3: continue
        out.append({"pts": [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], "x0": x0, "x1": x1, "y0": y0, "y1": y1, "rect": True})
    return out


def cap_height(page, known):
    """cap height (pt) of the plot-number font: from the first read's label runs whose text is a known plot number"""
    lib0 = [(c, tuple(b), a) for c, b, a in json.load(open(r"D:\code\FMB_to_GeoJSON\glyphlib.json"))]
    glyphs, runs, big, _ = F.build_glyphs(page, cls="upper"); hs = []
    for res in F.recognise(glyphs, runs, lib0):
        if res["is_label"] and (not known or res["text"] in known):
            for i in res["idx"]:
                P = [p for pt in glyphs[i]["parts"] for p in pt["pts"]]; hs.append(max(p[1] for p in P) - min(p[1] for p in P))
    return float(np.median(hs)) if len(hs) >= 2 else None


def read(pdf, scale, known=()):
    page = pdfplumber.open(pdf).pages[0]; H = float(page.height); M = F.PT_MM * scale / 1000.0
    hL = cap_height(page, set(known)) or cap_height(page, set())
    if hL is None: return [], None
    C = curves(page); dots = [((c["x0"] + c["x1"]) / 2, (c["y0"] + c["y1"]) / 2) for c in C if 0.5 <= max(c["x1"] - c["x0"], c["y1"] - c["y0"]) <= 1.6]
    G = []
    for i, o in enumerate(C):
        h = o["y1"] - o["y0"]; isJ = False
        if 1.2 * hL <= h <= 1.36 * hL and (o["x1"] - o["x0"]) / h < 0.34:            # capital J: cap height + descender
            bmJ, arJ = F._bitmap([o]); isJ = ham(bmJ, NARROW[1][1]) <= 12
            if isJ: G.append({"ch": "J", "d": ham(bmJ, NARROW[1][1]), "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y0"] + hL, "cy": o["y0"] + hL / 2})
            continue
        if abs(h - hL) > 0.045 * hL: continue
        poly = shapely.Polygon(o["pts"]).buffer(0)
        # true counters: small outlines lying wholly inside this character's outer ring
        parts = [o] + [q for j, q in enumerate(C) if j != i and (q["y1"] - q["y0"]) <= 0.93 * hL and max(q["x1"] - q["x0"], q["y1"] - q["y0"]) > 1.2
                       and q["x0"] > o["x0"] and q["x1"] < o["x1"] and q["y0"] > o["y0"] and q["y1"] < o["y1"]
                       and shapely.Polygon(q["pts"]).buffer(0).intersection(poly).area > 0.97 * max(shapely.Polygon(q["pts"]).buffer(0).area, 1e-9)]
        if o.get("rect"):
            if (o["x1"] - o["x0"]) / h < 0.22: G.append({"ch": "I", "d": 0, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "cy": (o["y0"] + o["y1"]) / 2})
            continue
        bm, ar = F._bitmap(parts); cand = {}
        for ch, b_, a_ in (NARROW if ar < 0.45 else LIB):
            if abs(math.log(ar / a_)) > 0.30: continue
            d = ham(bm, b_)
            if d < cand.get(ch, 999): cand[ch] = d
        if not cand: continue
        rk = sorted(cand.items(), key=lambda t: t[1]); best, bd = rk[0]; second = rk[1][1] if len(rk) > 1 else 999
        if bd > 36 or (bd > 22 and second - bd < 8): continue
        G.append({"ch": best, "d": bd, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "cy": (o["y0"] + o["y1"]) / 2})
    G.sort(key=lambda g: g["x0"]); par = list(range(len(G)))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a in range(len(G)):
        for b in range(a + 1, len(G)):
            if G[b]["x0"] - G[a]["x1"] > 0.6 * hL: break
            if abs(G[a]["cy"] - G[b]["cy"]) < 0.10 * hL and -0.12 * hL <= G[b]["x0"] - G[a]["x1"]: par[find(b)] = find(a)
    runs = {}
    for k in range(len(G)): runs.setdefault(find(k), []).append(G[k])
    out = []
    for r in runs.values():
        r.sort(key=lambda g: g["x0"]); x0, x1 = r[0]["x0"], r[-1]["x1"]; y0 = min(g["y0"] for g in r); y1 = max(g["y1"] for g in r)
        base = [(dx, dy) for dx, dy in dots if y1 - 0.35 * hL <= dy <= y1 + 0.2 * hL]
        inside = any(r[k]["x1"] - 0.05 * hL <= dx <= r[k + 1]["x0"] + 0.05 * hL for dx, dy in base for k in range(len(r) - 1))
        beside = any(x0 - 0.45 * hL <= dx <= x0 or x1 <= dx <= x1 + 0.45 * hL for dx, dy in base)
        txt = "".join(g["ch"] for g in r)
        dim = inside or (beside and txt.isdigit())
        strict = any(x0 - 0.75 * hL <= dx <= x1 + 0.75 * hL for dx, dy in base)
        out.append({"text": "".join(g["ch"] for g in r), "x": (x0 + x1) / 2 * M, "y": (H - (y0 + y1) / 2) * M, "worst": max(g["d"] for g in r), "isdim": dim, "strict": strict})
    return out, hL


if __name__ == "__main__":
    scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
    X = pd.read_csv(os.path.join(HERE, "out", "ledger.csv"), dtype=str).fillna(""); rows = []
    for (v, s), g in X.groupby(["village_code", "survey"], sort=False):
        S = gpd.read_file(os.path.join(HERE, f"out/fresh/{v}/{s}_parcels.geojson")); S["geometry"] = S.geometry.buffer(0)
        labs, hL = read(os.path.join(R, "FMB_Sketches", v, s + ".pdf"), scales.get((v, s), 500), known=[l for l in g.layer_plot if l])
        for L in labs:
            p = shapely.Point(L["x"], L["y"]); hit = [k for k, q in enumerate(S.geometry) if q.contains(p)]
            k = min(hit, key=lambda k: S.geometry.iloc[k].area) if hit else -1
            rows.append({"village_code": v, "survey": s, "text": L["text"], "worst": L["worst"], "isdim": L["isdim"], "strict": L["strict"], "x": round(L["x"], 2), "y": round(L["y"], 2), "pid": int(S.poly_id.iloc[k]) if k >= 0 else -1, "cap_pt": round(hL, 2)})
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "out", "labels3.csv"), index=False); print("labels:", len(rows))
