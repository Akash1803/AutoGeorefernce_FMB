"""Fourth read: the one-character-at-a-time reader of read3, for EVERY upright character size found on the sheet.
Portal sheets shrink the plot numbers on crowded layouts (390B: numbers 5.4 pt, dimensions 6.8 pt), so a reader tied to
one size leaves whole sheets unnumbered. A plot number has no decimal point; a dimension always has one."""
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
BLUE = (0.0, 0.0, 1.0)


def curves(page):
    out = []
    for c in page.curves:
        if c.get("non_stroking_color") != BLUE: continue
        pts = F.flatten(c["path"])
        if len(pts) < 3: continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        out.append({"pts": pts, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys)})
    for r in page.rects:                      # the letter I is a plain rectangle
        if r.get("non_stroking_color") != BLUE or not r.get("fill", True): continue
        x0, x1, y0, y1 = r["x0"], r["x1"], r["top"], r["bottom"]
        if x1 - x0 < 0.2 or y1 - y0 < 2.5: continue
        out.append({"pts": [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)], "x0": x0, "x1": x1, "y0": y0, "y1": y1, "rect": True})
    return out


def sizes(C):
    """upright character sizes on the sheet: heights shared by many outlines (rotated text spreads its heights)"""
    h = np.sort(np.array([c["y1"] - c["y0"] for c in C if 3.0 <= c["y1"] - c["y0"] <= 18.0 and not c.get("rect")])); out = []; i = 0
    while i < len(h):
        j = i
        while j + 1 < len(h) and h[j + 1] - h[i] <= 0.012 * h[i]: j += 1
        if j - i + 1 >= 3: out.append((float(np.median(h[i:j + 1])), j - i + 1))
        i = j + 1
    return out


def read(pdf, scale, force=()):
    page = pdfplumber.open(pdf).pages[0]; H = float(page.height); M = F.PT_MM * scale / 1000.0
    C = curves(page); S = sizes(C)
    if not S: return [], {}
    polys = {}
    def poly(i):
        if i not in polys: polys[i] = shapely.Polygon(C[i]["pts"]).buffer(0)
        return polys[i]
    order = sorted(range(len(C)), key=lambda i: C[i]["x0"]); xs0 = [C[i]["x0"] for i in order]
    def inside(i):
        o = C[i]; lo = bisect.bisect_left(xs0, o["x0"]); hi = bisect.bisect_right(xs0, o["x1"]); r = []
        for j in order[lo:hi]:
            q = C[j]
            if j != i and q["x0"] > o["x0"] and q["x1"] < o["x1"] and q["y0"] > o["y0"] and q["y1"] < o["y1"]: r.append(j)
        return r
    dots_all = [((c["x0"] + c["x1"]) / 2, (c["y0"] + c["y1"]) / 2, max(c["x1"] - c["x0"], c["y1"] - c["y0"])) for c in C if max(c["x1"] - c["x0"], c["y1"] - c["y0"]) <= 2.2]
    out = []; score = {}
    for hL, _n in S:
        G = []
        for i, o in enumerate(C):
            h = o["y1"] - o["y0"]
            if 1.2 * hL <= h <= 1.36 * hL and (o["x1"] - o["x0"]) / h < 0.34 and not o.get("rect"):            # capital J: cap height + descender
                bmJ, arJ = F._bitmap([o]); d = ham(bmJ, NARROW[1][1])
                if d <= 12: G.append({"ch": "J", "d": d, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y0"] + hL, "cy": o["y0"] + hL / 2})
                continue
            if abs(h - hL) > 0.03 * hL: continue
            if o.get("rect"):
                if (o["x1"] - o["x0"]) / h < 0.22: G.append({"ch": "I", "d": 0, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "cy": (o["y0"] + o["y1"]) / 2})
                continue
            parts = [o] + [C[j] for j in inside(i) if (C[j]["y1"] - C[j]["y0"]) <= 0.93 * hL and max(C[j]["x1"] - C[j]["x0"], C[j]["y1"] - C[j]["y0"]) > 0.16 * hL
                           and poly(j).intersection(poly(i)).area > 0.97 * max(poly(j).area, 1e-9)]
            bm, ar = F._bitmap(parts); cand = {}
            for ch, b_, a_ in (NARROW if ar < 0.45 else LIB):
                if abs(math.log(ar / a_)) > 0.30: continue
                d = ham(bm, b_)
                if d < cand.get(ch, 999): cand[ch] = d
            if not cand: continue
            rk = sorted(cand.items(), key=lambda t: t[1]); best, bd = rk[0]; second = rk[1][1] if len(rk) > 1 else 999
            if bd > 36 or (bd > 22 and second - bd < 8): continue
            G.append({"ch": best, "d": bd, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"], "cy": (o["y0"] + o["y1"]) / 2})
        if not G: continue
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
        dots = [(x, y) for x, y, m in dots_all if 0.07 * hL <= m <= 0.22 * hL]
        for r in runs.values():
            r.sort(key=lambda g: g["x0"]); x0, x1 = r[0]["x0"], r[-1]["x1"]; y0 = min(g["y0"] for g in r); y1 = max(g["y1"] for g in r)
            base = [(dx, dy) for dx, dy in dots if y1 - 0.35 * hL <= dy <= y1 + 0.2 * hL]
            ins = any(r[k]["x1"] - 0.05 * hL <= dx <= r[k + 1]["x0"] + 0.05 * hL for dx, dy in base for k in range(len(r) - 1))
            beside = any(x0 - 0.45 * hL <= dx <= x0 or x1 <= dx <= x1 + 0.45 * hL for dx, dy in base)
            txt = "".join(g["ch"] for g in r)
            out.append({"text": txt, "x": (x0 + x1) / 2 * M, "y": (H - (y0 + y1) / 2) * M, "px": (x0 + x1) / 2, "py": (y0 + y1) / 2, "worst": max(g["d"] for g in r), "isdim": ins or (beside and txt.isdigit()),
                        "strict": any(x0 - 0.75 * hL <= dx <= x1 + 0.75 * hL for dx, dy in base), "cap_pt": round(hL, 2)})
    # which sizes are plot-number sizes: the ones that spell clean numbers of two or more characters with no decimal point;
    # rotated dimension text lands in other bands as loose single digits and is left out
    for a in out:
        sc = score.setdefault(a["cap_pt"], {"multi": 0, "dim": 0, "all": 0}); sc["all"] += 1
        if a["isdim"]: sc["dim"] += 1
        elif len(a["text"]) >= 2 and a["worst"] <= 12 and not (a["text"].isdigit() and a["strict"]) and set(a["text"]) != {"0"}: sc["multi"] += 1
    top = max([v["multi"] for v in score.values()] or [0]); dimtop = max([v["dim"] for v in score.values()] or [0])
    label_sizes = {k for k, v in score.items() if v["multi"] >= 2 and v["multi"] >= 0.25 * top}
    dim_sizes = {k for k, v in score.items() if v["dim"] >= 3 and v["dim"] >= 0.5 * dimtop}
    if force: label_sizes |= {k for k in score if any(abs(k - f) <= 0.03 * f for f in force)}
    keep0 = [a for a in out if a["cap_pt"] in label_sizes or (a["cap_pt"] in dim_sizes and not label_sizes)]
    seen, keep = set(), []
    for a in sorted(keep0, key=lambda a: a["worst"]):           # the same characters can fall in two neighbouring bands: keep one copy
        k = (a["text"], round(a["px"]), round(a["py"]))
        if k not in seen: seen.add(k); keep.append(a)
    return keep, {"label_sizes": sorted(label_sizes), "dimension_sizes": sorted(dim_sizes), "bands": {k: v for k, v in score.items() if v["multi"] or v["dim"] >= 3}}


if __name__ == "__main__":
    import pandas as pd
    R = r"D:\Projects\Tambaram_Chengalpattu"; v, s = sys.argv[1], sys.argv[2]
    sc = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
    L, S = read(os.path.join(R, "FMB_Sketches", v, s + ".pdf"), sc.get((v, s), 500)); print(S)
    d = pd.DataFrame(L); lab = d[~d.isdim & ~(d.text.str.isdigit() & d.strict)]; print("labels", len(d), "| plot numbers", len(lab), "| by size", lab.cap_pt.value_counts().to_dict())
    print(sorted(lab.text.tolist(), key=lambda t: (len(t), t)))
