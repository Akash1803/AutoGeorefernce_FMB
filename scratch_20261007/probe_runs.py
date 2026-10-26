"""Probe: group black glyphs into runs, recognise them, and find underlines."""
import sys, math, collections
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import pdfplumber
import fmb_to_geojson as C
from probe_black import black_glyphs, templates

def glyphs_of(curves):
    curves = [c for c in curves if c["d"] > 1.2]
    curves.sort(key=lambda c: -c["d"])
    used, out = set(), []
    for i, o in enumerate(curves):
        if i in used: continue
        parts = [o]; used.add(i)
        for j, s in enumerate(curves):
            if j in used: continue
            if s["x0"] >= o["x0"] - .3 and s["x1"] <= o["x1"] + .3 and s["y0"] >= o["y0"] - .3 and s["y1"] <= o["y1"] + .3:
                parts.append(s); used.add(j)
        bm, ar = C._bitmap(parts)
        out.append({"parts": parts, "bmp": bm, "ar": ar, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"],
                    "h": o["y1"] - o["y0"], "cx": (o["x0"] + o["x1"]) / 2, "cy": (o["y0"] + o["y1"]) / 2})
    return out

def runs_of(gl):
    adj = collections.defaultdict(set)
    for i in range(len(gl)):
        for j in range(i + 1, len(gl)):
            a, b = gl[i], gl[j]
            hh = max(a["h"], b["h"])
            if abs(a["cy"] - b["cy"]) < 0.6 * hh and (abs(a["y1"] - b["y1"]) < 0.35 * hh or abs(a["cy"] - b["cy"]) < 0.25 * hh) \
               and max(a["x0"], b["x0"]) - min(a["x1"], b["x1"]) < 0.45 * hh:
                adj[i].add(j); adj[j].add(i)
    seen, runs = set(), []
    for i in range(len(gl)):
        if i in seen: continue
        st, c = [i], []
        while st:
            n = st.pop()
            if n in seen: continue
            seen.add(n); c.append(n); st += list(adj[n] - seen)
        runs.append(sorted(c, key=lambda k: gl[k]["x0"]))
    return runs

def ham(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def read(gl, run, T):
    txt, worst = "", 0
    for k in run:
        g = gl[k]; best, bd = "?", 999
        for ch, bm, ar in T:
            if abs(math.log(max(g["ar"], 1e-3) / ar)) > 0.45: continue
            d = ham(g["bmp"], bm)
            if d < bd: best, bd = ch, d
        txt += best if bd <= 60 else "?"; worst = max(worst, bd)
    return txt, worst

if __name__ == "__main__":
    s = sys.argv[1]
    page = pdfplumber.open(r"D:/Projects/Tambaram_Chengalpattu/FMB_Sketches/35_04_221/%s.pdf" % s).pages[0]
    H = float(page.height)
    gl = glyphs_of(black_glyphs(page)); T = templates()
    lines = [(l["x0"], l["x1"], l["top"]) for l in page.lines] + [(r["x0"], r["x1"], r["top"]) for r in page.rects if r["height"] < 2.5 and r["width"] > 3]
    for r in runs_of(gl):
        x0 = min(gl[k]["x0"] for k in r); x1 = max(gl[k]["x1"] for k in r); y1 = max(gl[k]["y1"] for k in r); y0 = min(gl[k]["y0"] for k in r)
        txt, worst = read(gl, r, T)
        ul = any(abs(lx0 - x0) < 6 and abs(lx1 - x1) < 6 and 0 < top - y1 < 8 for lx0, lx1, top in lines)
        print("%-14s worst=%3d  x=%7.1f y=%7.1f h=%4.1f  underlined=%s" % (txt, worst, (x0 + x1) / 2, (y0 + y1) / 2, y1 - y0, ul))
