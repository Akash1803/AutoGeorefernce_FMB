import sys, json, math, collections
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON"); import pdfplumber, fmb_to_geojson as F
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union
lib = [(c, tuple(b), x) for c, b, x in json.load(open(r"D:\code\FMB_to_GeoJSON\glyphlib.json"))]
W_, H_ = F.W_, F.H_

def black_runs(page):
    """glyph runs drawn in BLACK (neighbour survey numbers, title block), same grouping as build_glyphs."""
    curves = []
    for c in page.curves:
        if c.get("non_stroking_color") != (0.0, 0.0, 0.0): continue
        pts = F.flatten(c["path"])
        if len(pts) < 3: continue
        d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
        curves.append({"pts": pts, "d": d, "x0": min(p[0] for p in pts), "x1": max(p[0] for p in pts), "y0": min(p[1] for p in pts), "y1": max(p[1] for p in pts)})
    H = float(page.height); k = H / F.A0_H
    curves = [c for c in curves if 150 * k < H - (c["y0"] + c["y1"]) / 2 < H - 150 * k]      # drop title block / footer
    big = [c for c in curves if c["d"] >= 4.5]; small = [c for c in curves if c["d"] < 4.5]
    glyphs = []
    for o in big:
        parts = [o] + [s for s in small if s["d"] > 1.2 and s["x0"] >= o["x0"] - .3 and s["x1"] <= o["x1"] + .3 and s["y0"] >= o["y0"] - .3 and s["y1"] <= o["y1"] + .3]
        bm, ar = F._bitmap(parts); P = [p for pt in parts for p in pt["pts"]]
        glyphs.append({"parts": parts, "bmp": bm, "ar": ar, "d": o["d"], "cx": (o["x0"] + o["x1"]) / 2, "cy": (o["y0"] + o["y1"]) / 2, "y1": o["y1"]})
    med = sorted(g["d"] for g in glyphs)[len(glyphs) // 2] if glyphs else 8.0; sp = 1.1 * med
    adj = collections.defaultdict(set)
    for i in range(len(glyphs)):
        for j in range(i + 1, len(glyphs)):
            if math.dist((glyphs[i]["cx"], glyphs[i]["cy"]), (glyphs[j]["cx"], glyphs[j]["cy"])) <= sp: adj[i].add(j); adj[j].add(i)
    seen, runs = set(), []
    for i in range(len(glyphs)):
        if i in seen: continue
        st, c = [i], []
        while st:
            n = st.pop()
            if n in seen: continue
            seen.add(n); c.append(n); st += [m for m in adj[n] if m not in seen]
        runs.append(sorted(c, key=lambda k: glyphs[k]["cx"]))
    return glyphs, runs

def neighbour_labels(pdf):
    page = pdfplumber.open(pdf).pages[0]; H = float(page.height); M = F.PT_MM * 500 / 1000
    segs = F.extract_lines(page); black = [s for s in segs if s["layer"] in ("survey_boundary_main", "subdivision_line")]
    polys, _, _ = F.polygonize(black, 2.7 * H / F.A0_H)
    U = unary_union([Polygon(p["rings"][0]) for p in polys]).buffer(0)
    if U.geom_type == "MultiPolygon": U = max(U.geoms, key=lambda q: q.area)
    # underline candidates: short black horizontal lines anywhere (page.lines, y top-down)
    ul = [(l["x0"], l["x1"], l["top"]) for l in list(page.lines) + list(page.rects) if abs(l["top"] - l["bottom"]) < 1.5 and 4 < (l["x1"] - l["x0"]) < 60]
    glyphs, runs = black_runs(page); out = []
    for r in F.recognise(glyphs, runs, lib):
        txt = r["text"]
        if not txt or "?" in txt or not any(ch.isdigit() for ch in txt): continue
        gs = [glyphs[i] for i in r["idx"]]; x0 = min(g["cx"] for g in gs) - 3; x1 = max(g["cx"] for g in gs) + 3; ybot = max(g["y1"] for g in gs); gd = sum(g["d"] for g in gs) / len(gs)
        underlined = any(a <= x0 + 6 and b >= x1 - 6 and -0.2 * gd < t - ybot < 1.0 * gd for a, b, t in ul)
        cx = sum(g["cx"] for g in gs) / len(gs); cy = H - sum(g["cy"] for g in gs) / len(gs); pt = Point(cx, cy)
        inside = U.contains(pt); dist = U.exterior.distance(pt) * M
        c = U.centroid; bearing = (90 - math.degrees(math.atan2(cy - c.y, cx - c.x))) % 360
        out.append({"text": txt, "underlined": underlined, "inside": inside, "dist_m": round(dist, 1), "bearing": round(bearing), "pt": (cx, cy), "angle": r["angle"]})
    return out, U, polys, M, H

if __name__ == "__main__":
    S = r"C:\Users\FAI-AK~1\AppData\Local\Temp\claude\d--Akash-Python\be3d1523-acf3-4fef-9ca8-783aeaecb2d3\scratchpad"
    B = r"D:\Projects\Tambaram_Chengalpattu\FMB_Sketches"
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    for sheet in sys.argv[1:]:
        labs, U, polys, M, H = neighbour_labels(B + "\\" + sheet + ".pdf")
        nb = [l for l in labs if not l["inside"] and l["dist_m"] < 60]
        def compass(b): return ["N", "NE", "E", "SE", "S", "SW", "W", "NW"][int(((b + 22.5) % 360) // 45)]
        print(sheet, "| NEIGHBOURS (black, outside, <60 m):", [(l["text"], compass(l["bearing"]), l["dist_m"], "underlined" if l["underlined"] else "no-underline") for l in sorted(nb, key=lambda l: l["bearing"])])
        print("     other black runs seen:", [(l["text"], "in" if l["inside"] else "out", "ul" if l["underlined"] else "-") for l in labs if l not in nb][:14])
        fig, ax = plt.subplots(figsize=(8, 8))
        for p in polys:
            xs, ys = zip(*p["rings"][0]); ax.plot(xs, ys, color="k", lw=0.6)
        for l in labs:
            col = "red" if l in nb else ("grey" if not l["underlined"] else "orange")
            ax.text(l["pt"][0], l["pt"][1], l["text"], color=col, fontsize=9, ha="center", weight="bold" if l in nb else "normal")
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title("%s - red = neighbour survey numbers (underlined, outside)" % sheet)
        import os; fig.savefig(os.path.join(S, "nb_" + sheet.split("\\")[-1] + ".png"), dpi=110, bbox_inches="tight"); plt.close(fig)
