r"""
scale_from_dimensions.py
========================
Infer a sketch's true scale from its own edge dimensions.

Every boundary edge on an FMB sketch is annotated with its ground length
(e.g. 16.9). We read those annotations with the glyph library, pair each with
the nearest drawn segment, and compare annotated metres with drawn length at
the nominal scale. The median ratio, snapped to a standard scale, is the real
scale of the sheet. Independent of any external area figure.

    python scale_from_dimensions.py sketch.pdf [--nominal 500]
"""
import sys, os, json, math, statistics, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdfplumber
import fmb_to_geojson as F

STANDARD = [250, 400, 500, 800, 1000, 1250, 2000, 2500, 4000, 5000]


def dimension_runs(page, lib):
    """Yield (value_m, cx, cy, gd) for every recognised run that has a decimal point."""
    dots = []
    for c in page.curves:
        if c.get('non_stroking_color') != (0.0, 0.0, 1.0):
            continue
        pts = F.flatten(c['path'])
        if len(pts) < 3:
            continue
        xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if 0.5 <= max(w, h) <= 1.6 and h > 0 and 0.6 <= w / h <= 1.7:
            dots.append(((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2))
    bands = [b for b in F.size_groups(page) if b[1] >= 4.5]
    if not bands:
        return []
    band = (min(b[0] for b in bands), max(b[1] for b in bands))
    glyphs, runs, _, _ = F.build_glyphs(page, size_range=band)
    # Each dimension is split by its decimal point into a left run and a right
    # run (the dot is too small to be a glyph, and the gap breaks the run), so
    # stitch: for every dot find the run ending just left of it and the run
    # starting just right of it on the same baseline.
    R = []
    for res in F.recognise(glyphs, runs, lib):
        txt = res['text']
        if not txt or not txt.isdigit():
            continue
        gs = [glyphs[i] for i in res['idx']]
        R.append({'txt': txt, 'x0': gs[0]['cx'], 'x1': gs[-1]['cx'],
                  'cy': sum(g['cy'] for g in gs) / len(gs), 'gd': sum(g['d'] for g in gs) / len(gs), 'ang': res['angle']})
    out = []
    for dx, dy in dots:
        # work in each run's own frame (dimension text runs along its edge, so it
        # may be rotated): a = along-text coordinate, b = across-text offset
        def frame(r):
            th = math.radians(r['ang']); ux, uy = math.cos(th), math.sin(th)
            a0 = r['x0'] * ux + r['cy'] * uy; a1 = r['x1'] * ux + r['cy'] * uy
            ad = dx * ux + dy * uy; bd = -dx * uy + dy * ux; br = -r['x0'] * uy + r['cy'] * ux
            return a0, a1, ad, abs(bd - br)
        L, Rt = [], []
        for r in R:
            a0, a1, ad, off = frame(r)
            if off > 0.75 * r['gd']:
                continue
            if 0.2 * r['gd'] <= ad - max(a0, a1) <= 1.3 * r['gd']:
                L.append(r)
            elif 0.2 * r['gd'] <= min(a0, a1) - ad <= 1.3 * r['gd']:
                Rt.append(r)
        if len(L) != 1 or len(Rt) != 1:
            continue
        l, r = L[0], Rt[0]
        if abs(((l['ang'] - r['ang']) + 90) % 180 - 90) > 8:
            continue
        val = float(l['txt'] + '.' + r['txt'])
        out.append((val, (l['x0'] + r['x1']) / 2, (l['cy'] + r['cy']) / 2, l['gd'], l['ang'] if len(l['txt']) > 1 else r['ang']))
    return out


def estimate(pdf, lib, nominal=500):
    page = pdfplumber.open(pdf).pages[0]
    H = float(page.height); M = F.PT_MM * nominal / 1000.0
    segs = [s for s in F.extract_lines(page) if s['layer'] in ('survey_boundary_main', 'subdivision_line')]
    ratios = []
    for val, cx, cy, gd, ang in dimension_runs(page, lib):
        p = (cx, H - cy)
        best, bd = None, 1e9
        for s in segs:
            (ax, ay), (bx, by) = s['p0'], s['p1']
            L = math.hypot(bx - ax, by - ay)
            if L < 1e-6:
                continue
            seg_ang = math.degrees(math.atan2(by - ay, bx - ax)) % 180
            da = min(abs(seg_ang - ang % 180), 180 - abs(seg_ang - ang % 180))
            if da > 12:
                continue                                  # dimension text runs along its edge
            t = max(0, min(1, ((p[0] - ax) * (bx - ax) + (p[1] - ay) * (by - ay)) / (L * L)))
            d = math.hypot(p[0] - (ax + t * (bx - ax)), p[1] - (ay + t * (by - ay)))
            if d < bd:
                best, bd = L, d
        if best is None or bd > 2.5 * gd:
            continue
        drawn_m = best * M
        if drawn_m > 0.5:
            ratios.append(val / drawn_m)
    if len(ratios) < 3:
        return {'scale': nominal, 'confidence': 'low', 'n': len(ratios), 'median_ratio': None}
    r = statistics.median(ratios)
    scale = min(STANDARD, key=lambda s: abs(math.log(s / (nominal * r))))
    agree = sum(1 for x in ratios if 0.85 < x / r < 1.15) / len(ratios)
    return {'scale': scale, 'confidence': 'high' if agree >= 0.6 and len(ratios) >= 5 else 'medium',
            'n': len(ratios), 'median_ratio': round(r, 3), 'agree': round(agree, 2)}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('pdf'); ap.add_argument('--nominal', type=int, default=500); a = ap.parse_args()
    lib = [(c, tuple(b), x) for c, b, x in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'glyphlib.json')))]
    print(estimate(a.pdf, lib, a.nominal))
