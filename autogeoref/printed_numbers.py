"""Read the neighbour survey numbers printed around an FMB portal sheet, straight from the PDF.

Portal sheets have no text layer: every character is a filled outline. Plot numbers and edge
dimensions are blue (the converter reads those); the neighbouring survey numbers ("266", "68/A")
are black, horizontal, outside the parcel. They are drawn in DejaVu Sans, so each black glyph is
matched against DejaVu Sans bitmaps rendered here, on the converter's 14 x 20 grid.

A number is trusted as a neighbour only when the two sheets agree: sheet A prints B and sheet B
prints A (`agreed`). On Paranur, 2026-10-07, all 31 sheets were read in 33 s; the one-sided
pairs (264/71B, 70A/72, ...) are exactly the ones a reviewer has to settle.
"""
import json
import math
import os
import re
from dataclasses import asdict, dataclass

import numpy as np

from . import neighbours, paths

GRID_W, GRID_H = 14, 20             # the converter's glyph bitmap (fmb_to_geojson.W_, H_)
CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/-.()"
MAX_HAMMING = 40                    # worst glyph distance for a run to count as read
HEADER_SHARE, FOOTER_SHARE = 0.04, 0.972      # title block and footer bands, as page-height shares
PT_MM = 25.4 / 72.0
_TEMPLATES = None


@dataclass
class Label:
    text: str          # as printed, e.g. "68/A"
    number: str        # survey id it names, e.g. "68A" (neighbours.normalise)
    x: float           # sheet metres, the converter's frame (origin bottom-left)
    y: float
    worst: int         # worst glyph Hamming distance in the run


def _converter():
    import sys
    if str(paths.CODE) not in sys.path:
        sys.path.insert(0, str(paths.CODE))
    import fmb_to_geojson
    return fmb_to_geojson


def _flatten(path, n=8):
    return _converter().flatten(path, n)


def _bitmap(parts):
    return _converter()._bitmap(parts)


def templates():
    """[(char, bits, aspect)] for DejaVu Sans, rows ordered like the converter's bitmaps (y up)."""
    global _TEMPLATES
    if _TEMPLATES is None:
        import matplotlib
        from PIL import Image, ImageDraw, ImageFont
        font = ImageFont.truetype(os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts",
                                               "ttf", "DejaVuSans.ttf"), 200)
        out = []
        for ch in CHARS:
            img = Image.new("L", (400, 400), 0)
            ImageDraw.Draw(img).text((50, 50), ch, fill=255, font=font)
            a = np.array(img) > 127
            ys, xs = np.nonzero(a)
            if not len(xs):
                continue
            a = a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
            h, w = a.shape
            bits = []
            for r in range(GRID_H):
                y = min(int((GRID_H - 0.5 - r) * h / GRID_H), h - 1)       # row 0 = top of the glyph
                for c in range(GRID_W):
                    bits.append(1 if a[y, min(int((c + 0.5) * w / GRID_W), w - 1)] else 0)
            out.append((ch, tuple(bits), w / h))
        _TEMPLATES = out
    return _TEMPLATES


def _black_curves(page):
    out = []
    for c in page.curves:
        if c.get("non_stroking_color") != (0.0, 0.0, 0.0) or not c.get("fill"):
            continue
        pts = _flatten(c["path"])
        if len(pts) < 3:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
        if d > 1.2:
            out.append({"pts": pts, "d": d, "x0": min(xs), "x1": max(xs), "y0": min(ys), "y1": max(ys)})
    return out


def _glyphs(curves):
    """One glyph per outer outline, with the counters (holes of 0, 4, 6, 8, 9, A ...) it contains."""
    curves = sorted(curves, key=lambda c: -c["d"])
    used, out = set(), []
    for i, o in enumerate(curves):
        if i in used:
            continue
        parts = [o]
        used.add(i)
        for j, s in enumerate(curves):
            if j not in used and s["x0"] >= o["x0"] - .3 and s["x1"] <= o["x1"] + .3 \
                    and s["y0"] >= o["y0"] - .3 and s["y1"] <= o["y1"] + .3:
                parts.append(s)
                used.add(j)
        bm, ar = _bitmap(parts)
        out.append({"bmp": bm, "ar": ar, "x0": o["x0"], "x1": o["x1"], "y0": o["y0"], "y1": o["y1"],
                    "h": o["y1"] - o["y0"], "cy": (o["y0"] + o["y1"]) / 2})
    return out


def _runs(gl):
    """Glyphs on one text line, close enough to be one word."""
    adj = {i: set() for i in range(len(gl))}
    for i in range(len(gl)):
        for j in range(i + 1, len(gl)):
            a, b = gl[i], gl[j]
            hh = max(a["h"], b["h"])
            same_line = abs(a["cy"] - b["cy"]) < 0.6 * hh and (abs(a["y1"] - b["y1"]) < 0.35 * hh
                                                               or abs(a["cy"] - b["cy"]) < 0.25 * hh)
            if same_line and max(a["x0"], b["x0"]) - min(a["x1"], b["x1"]) < 0.45 * hh:
                adj[i].add(j)
                adj[j].add(i)
    seen, runs = set(), []
    for i in range(len(gl)):
        if i in seen:
            continue
        stack, comp = [i], []
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.append(n)
            stack += list(adj[n] - seen)
        runs.append(sorted(comp, key=lambda k: gl[k]["x0"]))
    return runs


def _read(gl, run):
    T = templates()
    txt, worst = "", 0
    for k in run:
        g = gl[k]
        best, bd = "?", 10 ** 6
        for ch, bm, ar in T:
            if abs(math.log(max(g["ar"], 1e-3) / ar)) > 0.45:
                continue
            d = sum(1 for x, y in zip(g["bmp"], bm) if x != y)
            if d < bd:
                best, bd = ch, d
        txt += best
        worst = max(worst, bd)
    return txt, worst


def read_labels(pdf_path, scale):
    """Every black survey-number label on the sheet, in sheet metres. Unreadable runs are dropped."""
    import pdfplumber
    page = pdfplumber.open(str(pdf_path)).pages[0]
    H = float(page.height)
    M = PT_MM * scale / 1000.0
    gl = _glyphs(_black_curves(page))
    out = []
    for run in _runs(gl):
        y0 = min(gl[k]["y0"] for k in run)
        y1 = max(gl[k]["y1"] for k in run)
        cy = (y0 + y1) / 2
        if cy < HEADER_SHARE * H or cy > FOOTER_SHARE * H:
            continue
        txt, worst = _read(gl, run)
        if worst > MAX_HAMMING or not re.fullmatch(r"\d+[A-Z]?(/[A-Z0-9]+)*", txt):
            continue
        number = neighbours.normalise(txt)
        if not number:
            continue
        cx = (min(gl[k]["x0"] for k in run) + max(gl[k]["x1"] for k in run)) / 2
        out.append(Label(txt, number, round(cx * M, 3), round((H - cy) * M, 3), int(worst)))
    return out


def sheet_scale(village, survey):
    f = paths.sheet_path(village, survey)
    props = json.loads(f.read_text(encoding="utf-8"))["features"][0]["properties"]
    return int(str(props.get("sheet_scale") or "1:500").split(":")[-1])


def village_table(village, cache_path=None):
    """{survey: [Label]} for every converted sheet of the village that has its PDF.

    Cached as JSON (keyed by each PDF's size and mtime) because reading 31 sheets takes ~30 s."""
    cache = {}
    if cache_path and os.path.exists(cache_path):
        cache = json.loads(open(cache_path, encoding="utf-8").read())
    table, changed = {}, False
    for s in paths.surveys_with_sheets(village):
        pdf = paths.sketch_dir(village) / ("%s.pdf" % s)
        if not pdf.exists():
            continue
        st = pdf.stat()
        key = "%d:%d" % (st.st_size, int(st.st_mtime))
        hit = cache.get(s)
        if hit and hit.get("key") == key:
            table[s] = [Label(**d) for d in hit["labels"]]
            continue
        table[s] = read_labels(pdf, sheet_scale(village, s))
        cache[s] = {"key": key, "labels": [asdict(L) for L in table[s]]}
        changed = True
    if cache_path and changed:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, indent=1)
    return table


def printed_sets(table):
    """{survey: set of survey ids it prints}, its own number left out."""
    return {s: {L.number for L in labels if L.number != str(s).upper()} for s, labels in table.items()}


def agreed(table, a, b):
    """True when both sheets are read and each prints the other."""
    p = printed_sets(table)
    a, b = str(a).upper(), str(b).upper()
    return a in p and b in p and b in p[a] and a in p[b]


def agreed_sets(table):
    """{survey: set of surveys it agrees with both ways}."""
    p = printed_sets(table)
    return {s: {n for n in p[s] if n in p and s.upper() in p[n]} for s in p}
