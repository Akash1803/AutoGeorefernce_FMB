"""Fifth read of the printed plot numbers: every blue character outline is classified once, by its own shape, whatever
its size; characters join into a number when they stand on one baseline, side by side, at the same size.
Why: portal sheets shrink the plot numbers on crowded layouts (390B: numbers 5.9 pt, dimensions 6.8 pt; 48: 4.5 pt), and a
reader tied to one size leaves such sheets unnumbered. A plot number has no decimal point; a dimension always has one
(the point is a small square, filed under rects, or a small round outline)."""
import sys, os, json, math, bisect, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, shapely, pdfplumber
import fmb_to_geojson as F
HERE = os.path.dirname(os.path.abspath(__file__))
FONT = json.load(open(os.path.join(HERE, "out", "font_bitmaps.json")))
BASE = [(c, tuple(b), a) for c, b, a in json.load(open(r"D:\code\FMB_to_GeoJSON\glyphlib.json")) if c not in "~()"]
LIB = BASE + [(c, tuple(FONT[c][0]), FONT[c][1]) for c in "0123456789ABCDEFGHKLMNPRSTUVWXYZ"]
NARROW = [(c, tuple(FONT[c][0]), FONT[c][1]) for c in "IJ"]
ham = lambda a, b: sum(1 for x, y in zip(a, b) if x != y)
BLUE = (0.0, 0.0, 1.0); JRATIO = 1.28


def outlines(page):
    C, dots = [], []
    for c in page.curves:
        if c.get("non_stroking_color") != BLUE: continue
        pts = F.flatten(c["path"])
        if len(pts) < 3: continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; o = {"pts": pts, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys)}
        m = max(o["x1"] - o["x0"], o["y1"] - o["y0"])
        if m <= 2.2: dots.append(((o["x0"] + o["x1"]) / 2, (o["y0"] + o["y1"]) / 2, m))
        C.append(o)
    for r in page.rects:
        if r.get("non_stroking_color") != BLUE or not r.get("fill", True): continue
        x0, x1, y0, y1 = r["x0"], r["x1"], r["top"], r["bottom"]; m = max(x1 - x0, y1 - y0)
        if m <= 2.2 and min(x1 - x0, y1 - y0) >= 0.3: dots.append(((x0 + x1) / 2, (y0 + y1) / 2, m)); continue
        if x1 - x0 < 0.2 or y1 - y0 < 2.5: continue
        C.append({"pts": [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], "x0": x0, "x1": x1, "y0": y0, "y1": y1, "rect": True})
    return C, dots


def read(pdf, scale):
    page = pdfplumber.open(pdf).pages[0]; H = float(page.height); M = F.PT_MM * scale / 1000.0
    C, dots = outlines(page)
    if not C: return [], {}
    polys = {}
    def poly(i):
        if i not in polys: polys[i] = shapely.Polygon(C[i]["pts"]).buffer(0)
        return polys[i]
    order = sorted(range(len(C)), key=lambda i: C[i]["x0"]); xs0 = [C[i]["x0"] for i in order]
    counter_of = {}
    cand = [i for i, o in enumerate(C) if 3.0 <= o["y1"] - o["y0"] <= 24.0]
    inner = {}
    for i in cand:
        o = C[i]; h = o["y1"] - o["y0"]; lo = bisect.bisect_left(xs0, o["x0"]); hi = bisect.bisect_right(xs0, o["x1"]); r = []
        if o.get("rect"): inner[i] = []; continue
        for j in order[lo:hi]:
            q = C[j]
            if j != i and q["x0"] > o["x0"] and q["x1"] < o["x1"] and q["y0"] > o["y0"] and q["y1"] < o["y1"] and (q["y1"] - q["y0"]) <= 0.93 * h and max(q["x1"] - q["x0"], q["y1"] - q["y0"]) > 0.16 * h \
                    and poly(j).intersection(poly(i)).area > 0.97 * max(poly(j).area, 1e-9):
                r.append(j); counter_of[j] = i
        inner[i] = r
    G = []
    for i in cand:
        if i in counter_of: continue                       # the hole of a 0, 6, 8, A ... is not a character of its own
        o = C[i]; h = o["y1"] - o["y0"]; w = o["x1"] - o["x0"]
        if o.get("rect"):
            if w / h < 0.22: G.append({"ch": "I", "d": 0, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "h": h})
            continue
        bm, ar = F._bitmap([o] + [C[j] for j in inner[i]]); c = {}
        for ch, b_, a_ in (NARROW if ar < 0.45 else LIB):
            if abs(math.log(ar / a_)) > 0.30: continue
            d = ham(bm, b_)
            if d < c.get(ch, 999): c[ch] = d
        if not c: continue
        rk = sorted(c.items(), key=lambda t: t[1]); best, bd = rk[0]; second = rk[1][1] if len(rk) > 1 else 999
        if bd > 36 or (bd > 22 and second - bd < 8): continue
        if best == "J":
            if bd > 12: continue
            he = h / JRATIO; G.append({"ch": "J", "d": bd, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y0"] + he, "h": he})
        elif best == "I": continue                         # a drawn I is a rectangle; a curved narrow outline is something else
        else: G.append({"ch": best, "d": bd, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "h": h})
    for g in G: g["cy"] = (g["y0"] + g["y1"]) / 2
    G.sort(key=lambda g: g["x0"]); par = list(range(len(G)))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for a in range(len(G)):
        ha = G[a]["h"]
        for b in range(a + 1, len(G)):
            if G[b]["x0"] - G[a]["x1"] > 0.6 * ha: break
            hb = G[b]["h"]; hm = (ha + hb) / 2
            if abs(ha - hb) <= 0.08 * hm and abs(G[a]["cy"] - G[b]["cy"]) < 0.10 * hm and abs(G[a]["y1"] - G[b]["y1"]) < 0.10 * hm and -0.12 * hm <= G[b]["x0"] - G[a]["x1"]: par[find(b)] = find(a)
    runs = {}
    for k in range(len(G)): runs.setdefault(find(k), []).append(G[k])
    out = []
    for r in runs.values():
        r.sort(key=lambda g: g["x0"]); x0, x1 = r[0]["x0"], r[-1]["x1"]; y0 = min(g["y0"] for g in r); y1 = max(g["y1"] for g in r); hL = float(np.median([g["h"] for g in r]))
        base = [(dx, dy) for dx, dy, m in dots if 0.07 * hL <= m <= 0.24 * hL and y1 - 0.35 * hL <= dy <= y1 + 0.2 * hL]
        ins = any(r[k]["x1"] - 0.05 * hL <= dx <= r[k + 1]["x0"] + 0.05 * hL for dx, dy in base for k in range(len(r) - 1))
        beside = any(x0 - 0.45 * hL <= dx <= x0 or x1 <= dx <= x1 + 0.45 * hL for dx, dy in base)
        txt = "".join(g["ch"] for g in r)
        out.append({"text": txt, "x": (x0 + x1) / 2 * M, "y": (H - (y0 + y1) / 2) * M, "px": (x0 + x1) / 2, "py": (y0 + y1) / 2, "worst": max(g["d"] for g in r), "isdim": ins or (beside and txt.isdigit()),
                    "strict": any(x0 - 0.75 * hL <= dx <= x1 + 0.75 * hL for dx, dy in base), "cap_pt": round(hL, 2)})
    # the sizes plot numbers are printed in on this sheet: taken from the clean numbers of two or more characters
    good = sorted(a["cap_pt"] for a in out if not a["isdim"] and len(a["text"]) >= 2 and a["worst"] <= 12 and not (a["text"].isdigit() and a["strict"]) and set(a["text"]) != {"0"})
    dim = sorted(a["cap_pt"] for a in out if a["isdim"])
    def peaks(v, least):
        P = []; i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[j + 1] - v[i] <= 0.06 * v[i]: j += 1
            P.append((float(np.median(v[i:j + 1])), j - i + 1)); i = j + 1
        top = max([n for _h, n in P] or [0]); return [(h, n) for h, n in P if n >= least and n >= 0.2 * top]
    LS = peaks(good, 2); DS = peaks(dim, 3)
    sizes = [h for h, _n in LS] or [h for h, _n in DS]            # a sheet of single-digit plots only: they are printed at the size of the dimensions
    for a in out: a["label_size"] = any(abs(a["cap_pt"] - h) <= 0.06 * h for h in sizes)
    if not LS:                                   # no longer number to take the size from: a plot number is printed 0.85 to 1.12 times the size of the dimensions
        for a in out:
            if not a["label_size"] and len(a["text"]) <= 2 and a["worst"] <= 8 and any(0.85 * h <= a["cap_pt"] <= 1.12 * h for h, _n in DS): a["label_size"] = True
    return out, {"label_sizes": [(round(h, 2), n) for h, n in LS], "dimension_sizes": [(round(h, 2), n) for h, n in DS]}


def plot_numbers(L):
    """the reads that can be plot numbers: right size, no decimal point, not a lone digit next to a point, clean shapes"""
    return [a for a in L if a["label_size"] and not a["isdim"] and not (a["text"].isdigit() and a["strict"]) and not (a["text"].isdigit() and len(a["text"]) > 3)
            and a["worst"] <= 22 and set(a["text"]) != {"0"}]


if __name__ == "__main__":
    import pandas as pd
    R = r"D:\Projects\Tambaram_Chengalpattu"; v, s = sys.argv[1], sys.argv[2]
    sc = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
    L, S = read(os.path.join(R, "FMB_Sketches", v, s + ".pdf"), sc.get((v, s), 500)); print(S)
    lab = plot_numbers(L); print("reads", len(L), "| plot numbers", len(lab))
    print(sorted([a["text"] for a in lab], key=lambda t: (len(t), t)))
