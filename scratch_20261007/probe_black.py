"""Probe: black glyph outlines on a portal FMB sheet (neighbour numbers), recognised against DejaVu Sans templates."""
import sys, math, collections, os
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import pdfplumber, numpy as np
import fmb_to_geojson as C
from PIL import Image, ImageDraw, ImageFont
import matplotlib
FONT = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf", "DejaVuSans.ttf")
CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/-.()"

def template(ch, font):
    img = Image.new("L", (400, 400), 0); d = ImageDraw.Draw(img); d.text((50, 50), ch, fill=255, font=font)
    a = np.array(img) > 127; ys, xs = np.nonzero(a)
    if not len(xs): return None
    a = a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]; h, w = a.shape
    bits = []
    for r in range(C.H_):
        for c in range(C.W_):
            y = int((C.H_ - 0.5 - r) * h / C.H_); x = int((c + 0.5) * w / C.W_)
            bits.append(1 if a[min(y, h - 1), min(x, w - 1)] else 0)
    return tuple(bits), w / h

def templates():
    f = ImageFont.truetype(FONT, 200)
    return [(ch,) + template(ch, f) for ch in CHARS if template(ch, f)]

def black_glyphs(page):
    curves = []
    for c in page.curves:
        if c.get("non_stroking_color") != (0.0, 0.0, 0.0) or not c.get("fill"): continue
        pts = C.flatten(c["path"])
        if len(pts) < 3: continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
        curves.append({"pts": pts, "d": d, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys)})
    return curves

if __name__ == "__main__":
    s = sys.argv[1]
    page = pdfplumber.open(r"D:/Projects/Tambaram_Chengalpattu/FMB_Sketches/35_04_221/%s.pdf" % s).pages[0]
    cv = black_glyphs(page)
    print("black filled curves:", len(cv), "sizes:", sorted(round(c["d"], 1) for c in cv)[:10], "...", sorted(round(c["d"], 1) for c in cv)[-10:])
    hs = collections.Counter(round(c["y1"] - c["y0"], 1) for c in cv); print("heights:", hs.most_common(10))
