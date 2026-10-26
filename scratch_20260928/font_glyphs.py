"""Exemplar bitmaps for A-Z / 0-9 from DejaVu Sans outlines, in the converter's own bitmap format (14x20, even-odd)."""
import sys, json, math; sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import fmb_to_geojson as F
from matplotlib import font_manager
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import DecomposingRecordingPen
def contours(font, ch, n=8):
    gs = font.getGlyphSet(); name = font.getBestCmap()[ord(ch)]; pen = DecomposingRecordingPen(gs); gs[name].draw(pen); out, cur, last = [], [], None
    for op, a in pen.value:
        if op == "moveTo": cur = [a[0]]; last = a[0]
        elif op == "lineTo": cur.append(a[0]); last = a[0]
        elif op == "qCurveTo":
            pts = list(a); on = []
            # TrueType: consecutive off-curve points have implied on-curve midpoints
            seq = [last] + pts
            k = 1
            while k < len(seq) - 1:
                c = seq[k]; nxt = seq[k + 1]
                end = nxt if k + 1 == len(seq) - 1 else ((c[0] + nxt[0]) / 2, (c[1] + nxt[1]) / 2)
                for t in range(1, n + 1):
                    u = t / n; cur.append(((1 - u) ** 2 * last[0] + 2 * (1 - u) * u * c[0] + u * u * end[0], (1 - u) ** 2 * last[1] + 2 * (1 - u) * u * c[1] + u * u * end[1]))
                last = end; k += 1
        elif op == "curveTo":
            c1, c2, end = a
            for t in range(1, n + 1):
                u = t / n; cur.append(((1 - u) ** 3 * last[0] + 3 * (1 - u) ** 2 * u * c1[0] + 3 * (1 - u) * u * u * c2[0] + u ** 3 * end[0], (1 - u) ** 3 * last[1] + 3 * (1 - u) ** 2 * u * c1[1] + 3 * (1 - u) * u * u * c2[1] + u ** 3 * end[1]))
            last = end
        elif op in ("closePath", "endPath"):
            if cur: cur.append(cur[0]); out.append({"pts": cur}); cur = []
    return out
if __name__ == "__main__":
    path = font_manager.findfont("DejaVu Sans"); font = TTFont(path); print(path)
    lib = json.load(open(r"D:\code\FMB_to_GeoJSON\glyphlib.json")); ham = lambda a, b: sum(1 for x, y in zip(a, b) if x != y)
    for ch in "0123456789ABC":
        bm, ar = F._bitmap(contours(font, ch)); ex = [(ham(bm, tuple(b)), round(a, 3)) for c, b, a in lib if c == ch]
        print(ch, "font aspect", round(ar, 3), "| hamming to library exemplars (of 280):", sorted(ex)[:4])
