#!/usr/bin/env python3
"""
fmb_to_geojson.py
=================

Convert a Tamil Nadu cadastral PDF -- either an FMB / survey sketch or a
CMDA / DTCP layout-approval plan -- into polygonised GeoJSON with plot-number
attributes.

The core idea: these PDFs are *vector* drawings. The parcel boundaries are real
line segments in the PDF, so we can recover true geometry (angles, lengths and
shapes intact) rather than tracing a raster. Plot numbers are recovered two
ways depending on the document:

  * FMB sketches have NO font layer -- the digits are drawn as vector glyph
    outlines. We read them by matching each glyph's shape against a library of
    known exemplars (see --glyphlib).
  * CMDA / DTCP layout plans DO have a real text layer, so plot numbers come
    straight from pdfplumber's extracted words.

Only two transforms are ever applied to the geometry, and both are conformal
(shape-preserving): a Y-flip (PDF y is top-down) and a single uniform scale to
ground metres using the sheet's stated scale. No rotation, snapping to a grid,
or simplification of the drawn shape. Node-snapping during polygonisation is
the one exception and is bounded by --snap (reported per run).

See fmb_to_geojson.md for the full method, the per-document-class notes, and
the known limitations.

USAGE
-----
    # FMB sketch (glyph-outline plot numbers)
    python fmb_to_geojson.py sketch.pdf --scale 500 --area 11853 \\
        --village "Koladi [70]" --taluk Poonamallee --survey 8 \\
        --glyphlib glyphlib.json -o out/

    # CMDA / DTCP layout plan (text-layer plot numbers), auto-detected
    python fmb_to_geojson.py layout.pdf --scale 800 --kind layout \\
        --clip 140 1160 1220 1780 -o out/

    # let the tool read scale/area/village from the title block first
    python fmb_to_geojson.py sketch.pdf --show-titleblock

Outputs (into --out, default ./out):
    <stem>_parcels.geojson   polygons with plot_no + area attributes
    <stem>_parcels.csv       the attribute table, sorted by plot_no
    <stem>_lines.geojson     every classified line segment
    <stem>_preview.png       labelled parcel map for visual QC

Dependencies: pdfplumber, matplotlib (preview only). No network needed.
"""

import argparse, json, math, csv, os, sys, collections, itertools

PT_MM = 25.4 / 72.0            # one PostScript point in millimetres
A0_H = 3370.0                  # reference sheet height; thresholds scale off this


# ----------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------
def flatten(path, n=8):
    """Flatten a pdfplumber path (m/l/v/h ops) to a list of points."""
    pts, cur, start = [], None, None
    for t in path:
        op = t[0]
        if op == 'm':
            cur = t[1]; start = cur; pts.append(cur)
        elif op == 'l':
            pts.append(t[1]); cur = t[1]
        elif op == 'v':                       # cubic bezier: (c2, end), c1=cur
            c1, c2, e = cur, t[1], t[2]
            for i in range(1, n + 1):
                s = i / n; m = 1 - s
                pts.append((m*m*m*cur[0] + 3*m*m*s*c1[0] + 3*m*s*s*c2[0] + s*s*s*e[0],
                            m*m*m*cur[1] + 3*m*m*s*c1[1] + 3*m*s*s*c2[1] + s*s*s*e[1]))
            cur = e
        elif op == 'h' and start:
            pts.append(start)
    return pts


def shoelace(r):
    A = 0.0
    for i in range(len(r)):
        x1, y1 = r[i]; x2, y2 = r[(i + 1) % len(r)]
        A += x1 * y2 - x2 * y1
    return A / 2.0


def proj(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return None
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2
    return t, math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy))


def seg_x(p1, p2, p3, p4):
    d1x, d1y = p2[0] - p1[0], p2[1] - p1[1]
    d2x, d2y = p4[0] - p3[0], p4[1] - p3[1]
    den = d1x * d2y - d1y * d2x
    if abs(den) < 1e-12:
        return None
    t = ((p3[0] - p1[0]) * d2y - (p3[1] - p1[1]) * d2x) / den
    u = ((p3[0] - p1[0]) * d1y - (p3[1] - p1[1]) * d1x) / den
    if 1e-9 < t < 1 - 1e-9 and 1e-9 < u < 1 - 1e-9:
        return (p1[0] + t * d1x, p1[1] + t * d1y)
    return None


# ----------------------------------------------------------------------------
# line extraction + classification
# ----------------------------------------------------------------------------
def is_furniture(ln, W, H):
    """True if the line is page border / title-block, by frame geometry.

    Position- and orientation-based, NOT length-based: an unsubdivided survey
    can have one straight boundary hundreds of points long, and some sheets
    overflow the page, so neither length nor edge-proximity alone is safe.
    Thresholds scale with sheet height so A0 and A4 prints behave the same.
    """
    (x0, y0) = ln['path'][0][1]; (x1, y1) = ln['path'][-1][1]
    col = tuple(round(c, 3) for c in (ln.get('stroking_color') or (0, 0, 0)))
    if col not in ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)):
        return True                                    # grey frame strokes
    k = H / A0_H
    if max(y0, y1) < 150 * k or min(y0, y1) > H - 150 * k:
        return True                                    # title/footer band
    L = math.hypot(x1 - x0, y1 - y0)
    if L > 500 * k:
        if abs(y1 - y0) < 1.0 and (min(y0, y1) < 80 * k or max(y0, y1) > H - 80 * k):
            return True
        if abs(x1 - x0) < 1.0 and (min(x0, x1) < 80 * k or max(x0, x1) > W - 80 * k):
            return True
    return False


def extract_lines(page, clip=None):
    """Return classified sketch line segments, y-flipped to origin-bottom-left.

    Layers, by FMB drawing convention:
      survey_boundary_main  heaviest solid black stroke on the sheet
      subdivision_line      lighter solid black -- interior parcel divisions
      tie_line              dash-dot-dot -- surveyor check measurements (NOT
                            boundaries; excluded from polygonisation)
      blue_feature_line     blue -- natural features (channels, flow arrows)
    On layout plans every line is uniform black, so all solid black lines
    become subdivision_line and the classifier simply doesn't populate the
    tie/blue layers.
    """
    W, H = float(page.width), float(page.height)

    def inclip(x, ytop):
        if clip is None:
            return True
        y = H - ytop
        return clip[0] <= x <= clip[1] and clip[2] <= y <= clip[3]

    # heaviest solid-black stroke width on this sheet == the main boundary
    lws = []
    for ln in page.lines:
        if is_furniture(ln, W, H):
            continue
        col = tuple(round(c, 3) for c in (ln.get('stroking_color') or (0, 0, 0)))
        dash = (ln.get('dash') or ([], 0))[0]
        if col == (0.0, 0.0, 0.0) and len(dash) <= 4:
            lws.append(round(ln.get('linewidth') or 1.0, 1))
    heavy = max(lws) if lws else 1.0

    out = []
    for ln in page.lines:
        if is_furniture(ln, W, H):
            continue
        (x0, yt0) = ln['path'][0][1]; (x1, yt1) = ln['path'][-1][1]
        if math.hypot(x1 - x0, yt1 - yt0) < 1e-6:
            continue
        if clip is not None and not (inclip(x0, yt0) and inclip(x1, yt1)):
            continue
        col = tuple(round(c, 3) for c in (ln.get('stroking_color') or (0, 0, 0)))
        lw = round(ln.get('linewidth') or 1.0, 1)
        dash = (ln.get('dash') or ([], 0))[0]
        dashed = bool(dash) and len(dash) > 1 and dash[1] != 0
        if len(dash) > 4:
            layer = 'tie_line'
        elif col == (0.0, 0.0, 1.0):
            layer = 'blue_feature_line'
        elif heavy > 1.0 and lw >= heavy:
            layer = 'survey_boundary_main'
        else:
            layer = 'subdivision_line'
        out.append({'p0': (x0, H - yt0), 'p1': (x1, H - yt1),
                    'layer': layer, 'lw': lw, 'dashed': dashed})
    return out


# ----------------------------------------------------------------------------
# annotation arrows and slivers (2026-09-21, after 47A plot 3)
# ----------------------------------------------------------------------------
HEAD_ANGLE_DEG = (4.0, 60.0)     # an arrowhead stroke leans back along its line by this much
SLIVER_M2 = 2.5                  # faces smaller than this are folded into a neighbour


def drop_arrows(segs, tol):
    """Remove offset arrows from the polygonisation input; returns (kept, dropped).

    FMB sheets draw offset measurements to reference stones as thin black lines with a
    two-stroke arrowhead at the far end, the foot just inside the parcel. Polygonised, the
    foot snaps onto a nearby line and the arrow cuts a sliver off a real plot (47A plot 3,
    2026-09-21). An arrow is a line with at least two short strokes attached at one end,
    each leaning back along the line by HEAD_ANGLE_DEG, on opposite sides of it.
    """
    head_max = 5.0 * tol
    short = [i for i, s in enumerate(segs) if math.dist(s['p0'], s['p1']) <= head_max]
    drop = set()
    for i, s in enumerate(segs):
        L = math.dist(s['p0'], s['p1'])
        if L <= 1.5 * head_max or i in drop:
            continue
        for tip, other in ((s['p0'], s['p1']), (s['p1'], s['p0'])):
            ux, uy = (other[0] - tip[0]) / L, (other[1] - tip[1]) / L
            heads, sides, far = [], [], []
            for j in short:
                if j == i or j in drop:
                    continue
                h = segs[j]
                for a, b in ((h['p0'], h['p1']), (h['p1'], h['p0'])):
                    if math.dist(a, tip) > tol:
                        continue
                    hl = math.dist(a, b)
                    if hl < 1e-6:
                        break
                    vx, vy = (b[0] - a[0]) / hl, (b[1] - a[1]) / hl
                    ang = math.degrees(math.acos(max(-1.0, min(1.0, ux * vx + uy * vy))))
                    if HEAD_ANGLE_DEG[0] <= ang <= HEAD_ANGLE_DEG[1]:
                        heads.append(j); sides.append(ux * vy - uy * vx > 0); far.append(b)
                    break
            if len(heads) >= 2 and len(set(sides)) == 2:
                drop.add(i); drop.update(heads)
                for j in short:                              # the bar closing the head, if drawn
                    if j in drop:
                        continue
                    h = segs[j]
                    if (any(math.dist(h['p0'], f) <= tol for f in far)
                            and any(math.dist(h['p1'], f) <= tol for f in far)):
                        drop.add(j)
                break
    return [s for i, s in enumerate(segs) if i not in drop], [segs[i] for i in sorted(drop)]


def merge_slivers(polys, min_area=SLIVER_M2):
    """Fold faces smaller than `min_area` (same units as `area`) into the neighbour that shares
    the longest edge with them; returns (polys, merged). A sliver with no neighbour stays."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    geoms = [Polygon(p['rings'][0], p['rings'][1:]).buffer(0) for p in polys]
    merged, changed = 0, True
    while changed:
        changed = False
        for i, g in enumerate(geoms):
            if g is None or g.is_empty or g.area >= min_area:
                continue
            best, bl = None, 0.0
            for j, h in enumerate(geoms):
                if j == i or h is None or h.is_empty:
                    continue
                shared = g.boundary.intersection(h.boundary).length
                if shared > bl:
                    best, bl = j, shared
            if best is None or bl < 1e-6:
                continue
            u = unary_union([geoms[best], g]).buffer(0)
            if u.geom_type != 'Polygon':
                continue
            geoms[best], geoms[i] = u, None
            merged += 1; changed = True
    out = []
    for g in geoms:
        if g is None or g.is_empty:
            continue
        rings = [[list(c) for c in g.exterior.coords]] + [[list(c) for c in r.coords] for r in g.interiors]
        out.append({'area': g.area, 'rings': rings})
    out.sort(key=lambda p: -p['area'])
    return out, merged


# ----------------------------------------------------------------------------
# polygonisation (planar face traversal)
# ----------------------------------------------------------------------------
def polygonize(segs, tol):
    """Turn boundary segments into closed polygons.

    Steps: split at crossings and T-junctions, snap nodes within `tol`, prune
    dangling (degree-1) edges, walk the planar half-edge graph into faces, drop
    the unbounded outer face, and re-attach nested islands (separate components
    enclosed by a face) as holes rather than double-counting them.
    Returns (polys, dangles, stats) with geometry in the segments' own units.
    """
    S = [(s['p0'], s['p1']) for s in segs]
    cuts = collections.defaultdict(set)
    for i, j in itertools.combinations(range(len(S)), 2):
        pt = seg_x(*S[i], *S[j])
        if pt:
            for k in (i, j):
                cuts[k].add(round(proj(pt, *S[k])[0], 9))
    for i, (a, b) in enumerate(S):
        for j, (c, d) in enumerate(S):
            if i == j:
                continue
            for p in (c, d):
                r = proj(p, a, b)
                if r and r[1] <= tol and 1e-6 < r[0] < 1 - 1e-6:
                    cuts[i].add(round(r[0], 9))

    pieces = []
    for i, (a, b) in enumerate(S):
        ts = sorted({0.0, 1.0} | cuts.get(i, set()))
        for t0, t1 in zip(ts, ts[1:]):
            p0 = (a[0] + t0 * (b[0] - a[0]), a[1] + t0 * (b[1] - a[1]))
            p1 = (a[0] + t1 * (b[0] - a[0]), a[1] + t1 * (b[1] - a[1]))
            if math.dist(p0, p1) > 1e-6:
                pieces.append((p0, p1))

    grid, reps, shift = {}, [], 0.0

    def node(p):
        nonlocal shift
        gx, gy = int(p[0] // tol), int(p[1] // tol)
        best, bd = None, tol
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for k in grid.get((gx + dx, gy + dy), ()):
                    dd = math.dist(p, reps[k])
                    if dd < bd:
                        best, bd = k, dd
        if best is not None:
            shift = max(shift, bd)
            return best
        reps.append(p); grid.setdefault((gx, gy), []).append(len(reps) - 1)
        return len(reps) - 1

    adj = collections.defaultdict(set)
    for p0, p1 in pieces:
        u, v = node(p0), node(p1)
        if u != v:
            adj[u].add(v); adj[v].add(u)

    pruned, dangles = 0, []
    while True:
        dg = [n for n in adj if len(adj[n]) == 1]
        if not dg:
            break
        for n in dg:
            for m in list(adj[n]):
                dangles.append((reps[n], reps[m])); adj[m].discard(n)
            del adj[n]; pruned += 1
    adj = {n: v for n, v in adj.items() if v}
    if not adj:
        return [], [], {'nodes': 0, 'shift': shift, 'pruned': pruned,
                        'outer': 0.0, 'components': 0}

    comp, cid = {}, 0
    for s in adj:
        if s in comp:
            continue
        st, cid = [s], cid + 1
        while st:
            n = st.pop()
            if n in comp:
                continue
            comp[n] = cid; st += [m for m in adj[n] if m not in comp]

    order = {n: sorted(adj[n], key=lambda m: math.atan2(reps[m][1] - reps[n][1],
                                                        reps[m][0] - reps[n][0]))
             for n in adj}
    seen, faces = set(), []
    for u in adj:
        for v in adj[u]:
            if (u, v) in seen:
                continue
            ring, cu, cv = [], u, v
            while (cu, cv) not in seen:
                seen.add((cu, cv)); ring.append(cu)
                nb = order[cv]; cu, cv = cv, nb[(nb.index(cu) - 1) % len(nb)]
            if len(ring) >= 3:
                faces.append(ring)

    signed = [(shoelace([reps[i] for i in r]), r) for r in faces]
    negs = sorted([(A, r) for A, r in signed if A < 0])
    outer = negs[0] if negs else (0.0, [])
    islands = negs[1:]
    interior = [[A, r] for A, r in signed if A >= 1.0]

    def inring(ring, p):
        ins = False; n = len(ring)
        for i in range(n):
            x1, y1 = reps[ring[i]]; x2, y2 = reps[ring[(i + 1) % n]]
            if (y1 > p[1]) != (y2 > p[1]):
                if x1 + (p[1] - y1) * (x2 - x1) / (y2 - y1) > p[0]:
                    ins = not ins
        return ins

    holes = collections.defaultdict(list)
    for A, r in islands:
        probe = reps[r[0]]; mine = comp[r[0]]
        cand = [k for k, (fa, fr) in enumerate(interior)
                if comp[fr[0]] != mine and inring(fr, probe)]
        if not cand:
            interior.append([abs(A), list(reversed(r))]); continue
        h = min(cand, key=lambda k: interior[k][0])
        holes[h].append((abs(A), r))

    polys = []
    for k, (A, r) in enumerate(interior):
        net = A - sum(ha for ha, _ in holes.get(k, []))
        rings = [[list(reps[i]) for i in r] + [list(reps[r[0]])]]
        for _, hr in holes.get(k, []):
            rings.append([list(reps[i]) for i in hr] + [list(reps[hr[0]])])
        polys.append({'area': net, 'rings': rings})
    polys.sort(key=lambda p: -p['area'])
    return polys, dangles, {'nodes': len(adj), 'shift': shift, 'pruned': pruned,
                            'outer': abs(outer[0]), 'components': cid}


# ----------------------------------------------------------------------------
# glyph-outline plot numbers (FMB sheets with no font layer)
# ----------------------------------------------------------------------------
W_, H_ = 14, 20

def _bitmap(parts):
    P = [p for pt in parts for p in pt['pts']]
    x0, x1 = min(p[0] for p in P), max(p[0] for p in P)
    y0, y1 = min(p[1] for p in P), max(p[1] for p in P)
    sx, sy = (x1 - x0) or 1e-6, (y1 - y0) or 1e-6
    bits = []
    for r in range(H_):
        yy = y0 + sy * (H_ - 0.5 - r) / H_
        xs = []
        for pt in parts:
            q = pt['pts']
            for i in range(len(q) - 1):
                (ax, ay), (bx, by) = q[i], q[i + 1]
                if (ay > yy) != (by > yy):
                    xs.append(ax + (yy - ay) * (bx - ax) / (by - ay))
        for c in range(W_):
            xx = x0 + sx * (c + 0.5) / W_
            bits.append(1 if sum(1 for v in xs if v < xx) % 2 else 0)
    return tuple(bits), sx / sy


def _kmeans2(vals):
    lo, hi = min(vals), max(vals)
    for _ in range(80):
        mid = (lo + hi) / 2
        a = [d for d in vals if d < mid]; b = [d for d in vals if d >= mid]
        if not a or not b:
            return None
        na, nb = sum(a) / len(a), sum(b) / len(b)
        if abs((na + nb) / 2 - mid) < 1e-4:
            break
        lo, hi = na, nb
    return (lo + hi) / 2, lo, hi


def build_glyphs(page, cls='upper', spf=1.1, size_range=None):
    """Group blue glyph outlines into character runs (plot-number text).

    Two text sizes coexist on a sheet: plot numbers and edge dimensions. Which
    size the plot numbers are is not constant (larger on Koladi sheets, smaller
    on some others), so `cls` selects 'upper' or 'lower' size class. A 1-D
    k-means on glyph diameter finds the split.
    """
    curves = []
    for c in page.curves:
        if c.get('non_stroking_color') != (0.0, 0.0, 1.0):
            continue
        pts = flatten(c['path'])
        if len(pts) < 3:
            continue
        d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
        curves.append({'pts': pts, 'd': d,
                       'x0': min(p[0] for p in pts), 'x1': max(p[0] for p in pts),
                       'y0': min(p[1] for p in pts), 'y1': max(p[1] for p in pts)})
    if not curves:
        return [], [], 0.0, []

    allc = [c['d'] for c in curves]
    if size_range is not None:
        # explicit size band (see size_groups): glyphs of this text size only
        big, hi_cut = size_range[0] - 1e-6, size_range[1] + 1e-6
    else:
        r1 = _kmeans2(allc) if len(allc) >= 8 else None
        c1 = r1[0] if r1 else 4.5
        text = [d for d in allc if d >= c1]
        r2 = _kmeans2(text) if len(text) >= 8 else None
        hi_cut = None
        if r2 and r2[2] / r2[1] > 1.15:
            big = c1 if cls == 'lower' else r2[0]
            hi_cut = r2[0] if cls == 'lower' else None
        else:
            big = c1

    outers = [c for c in curves if c['d'] >= big and (hi_cut is None or c['d'] < hi_cut)]
    small = [c for c in curves if c['d'] < big]
    used, glyphs = set(), []
    for o in outers:
        parts = [o]
        for i, s in enumerate(small):
            if i in used or s['d'] <= 1.2:
                continue
            if (s['x0'] >= o['x0'] - .3 and s['x1'] <= o['x1'] + .3 and
                    s['y0'] >= o['y0'] - .3 and s['y1'] <= o['y1'] + .3):
                parts.append(s); used.add(i)
        P = [p for pt in parts for p in pt['pts']]
        bm, ar = _bitmap(parts)
        glyphs.append({'parts': parts, 'bmp': bm, 'ar': ar, 'd': o['d'],
                       'cx': (min(p[0] for p in P) + max(p[0] for p in P)) / 2,
                       'cy': (min(p[1] for p in P) + max(p[1] for p in P)) / 2})

    med = sorted(g['d'] for g in glyphs)[len(glyphs) // 2] if glyphs else 8.0
    sp = spf * med
    adj = collections.defaultdict(set)
    for i in range(len(glyphs)):
        for j in range(i + 1, len(glyphs)):
            if math.dist((glyphs[i]['cx'], glyphs[i]['cy']),
                         (glyphs[j]['cx'], glyphs[j]['cy'])) <= sp:
                adj[i].add(j); adj[j].add(i)
    seen, runs = set(), []
    for i in range(len(glyphs)):
        if i in seen:
            continue
        st, c = [i], []
        while st:
            n = st.pop()
            if n in seen:
                continue
            seen.add(n); c.append(n); st += [m for m in adj[n] if m not in seen]
        runs.append(sorted(c, key=lambda k: glyphs[k]['cx']))
    return glyphs, runs, big, []


def size_groups(page, min_count=4):
    """Distinct text sizes on the sheet, as (lo, hi) glyph-diameter bands.

    Portal sheets carry up to three sizes (tiny marks, edge dimensions, plot
    numbers) that can differ by only ~6 %, so bands are found by splitting the
    sorted diameters wherever there is a gap larger than 4 % (min 0.25 pt).
    """
    ds = []
    for c in page.curves:
        if c.get('non_stroking_color') != (0.0, 0.0, 1.0):
            continue
        pts = flatten(c['path'])
        if len(pts) < 3:
            continue
        d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
        if d >= 2.0:
            ds.append(d)
    ds.sort()
    if not ds:
        return []
    groups, start, prev, n = [], ds[0], ds[0], 1
    for d in ds[1:]:
        if d - prev > max(0.25, 0.04 * prev):
            groups.append((start, prev, n)); start, n = d, 0
        prev = d; n += 1
    groups.append((start, prev, n))
    return [(round(lo, 2), round(hi, 2)) for lo, hi, n in groups if n >= min_count]


def recognise(glyphs, runs, lib):
    """Classify each run of glyphs into text using a nearest-bitmap match.

    A run is a plot-number label only if it is horizontal (dimensions run
    tilted along their edge), contains no bracket (bracketed values are tie
    distances), and every glyph matched an exemplar.
    """
    def ham(a, b):
        return sum(1 for x, y in zip(a, b) if x != y)
    out = []
    for r in runs:
        C = [(glyphs[i]['cx'], glyphs[i]['cy']) for i in r]
        if len(C) > 1:
            mx = sum(c[0] for c in C) / len(C); my = sum(c[1] for c in C) / len(C)
            sxx = sum((c[0] - mx) ** 2 for c in C); syy = sum((c[1] - my) ** 2 for c in C)
            sxy = sum((c[0] - mx) * (c[1] - my) for c in C)
            ang = math.degrees(0.5 * math.atan2(2 * sxy, sxx - syy))
        else:
            ang = 0.0
        chars, worst = [], 0
        for i in r:
            g = glyphs[i]
            if g['ar'] < 0.45:            # brackets: far narrower than any digit
                chars.append('('); continue
            best, bd = '?', 999
            for ch, bm, ar in lib:
                if abs(math.log(g['ar'] / ar)) > 0.35:
                    continue
                dd = ham(g['bmp'], bm)
                if dd < bd:
                    best, bd = ch, dd
            if best == '~' and bd <= 30:
                continue                   # inner counter of 0/6/8/9, not a character
            chars.append(best if bd <= 30 else '?')
            worst = max(worst, bd)
        txt = ''.join(chars)
        is_label = (abs(ang) < 5.0 and txt != '' and '?' not in txt
                    and '(' not in txt and ')' not in txt)
        out.append({'idx': r, 'text': txt, 'is_label': is_label,
                    'angle': round(ang, 2), 'worst': worst})
    return out


# ----------------------------------------------------------------------------
# spatial join + output
# ----------------------------------------------------------------------------
def _inring(ring, p):
    ins = False
    for i in range(len(ring) - 1):
        x1, y1 = ring[i]; x2, y2 = ring[i + 1]
        if (y1 > p[1]) != (y2 > p[1]):
            if x1 + (p[1] - y1) * (x2 - x1) / (y2 - y1) > p[0]:
                ins = not ins
    return ins


def _inpoly(rings, p):
    if not _inring(rings[0], p):
        return False
    return not any(_inring(h, p) for h in rings[1:])


def run(pdf_path, scale, out_dir, kind='auto', glyphlib=None, area=None,
        village=None, taluk=None, survey=None, clip=None, snap_pt=None,
        cls='upper'):
    import pdfplumber
    page = pdfplumber.open(pdf_path).pages[0]
    H = float(page.height)
    M = PT_MM * scale / 1000.0                       # ground metres per point
    tol = (snap_pt if snap_pt is not None else 2.7) * H / A0_H

    if kind == 'auto':
        kind = 'layout' if len(page.chars) > 2000 else 'sketch'

    segs = extract_lines(page, clip=clip)
    black = [s for s in segs if s['layer'] in ('survey_boundary_main', 'subdivision_line')]
    black, arrows = drop_arrows(black, tol)
    for a_ in arrows:
        a_['layer'] = 'annotation_arrow'            # kept in the lines file, out of the polygons
    n_arrows = len(arrows)
    polys, dangles, st = polygonize(black, tol)
    for p in polys:
        p['area'] *= M * M
        p['rings'] = [[[round(x * M, 4), round(y * M, 4)] for x, y in r] for r in p['rings']]
    polys, n_slivers = merge_slivers(polys, SLIVER_M2)
    st['arrows_dropped'], st['slivers_merged'] = n_arrows, n_slivers

    # --- plot numbers ---
    labels = []          # each: {text, pt(metres), kind}
    if kind == 'layout':
        for w in page.extract_words():
            t = w['text'].strip()
            cx = (w['x0'] + w['x1']) / 2; cy = H - (w['top'] + w['bottom']) / 2
            if clip and not (clip[0] <= cx <= clip[1] and clip[2] <= (cy) <= clip[3]):
                continue
            up = t.upper()
            if t.replace('.', '').isdigit() or up in ('PARK', 'PP1', 'PP2', 'OSR', 'EWS'):
                labels.append({'text': up, 'pt': (cx * M, cy * M), 'ambiguous': False})
    else:
        lib = None
        if glyphlib and os.path.exists(glyphlib):
            lib = [(c, tuple(b), a) for c, b, a in json.load(open(glyphlib))]
        if lib:
            # decimal points: tiny round blue blobs. Every edge dimension on a
            # portal sheet carries one (16.9, 9.2 ...); plot numbers never do.
            dots = []
            for c in page.curves:
                if c.get('non_stroking_color') != (0.0, 0.0, 1.0):
                    continue
                pts = flatten(c['path'])
                if len(pts) < 3:
                    continue
                d = max(math.dist(a, b) for a in pts[::3] for b in pts[::3])
                xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
                w, h = max(xs) - min(xs), max(ys) - min(ys)
                if 0.5 <= max(w, h) <= 1.6:
                    if h > 0 and 0.6 <= w / h <= 1.7:
                        dots.append(((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2))

            def _labels_for(c, band=None):
                out = []
                glyphs, runs, big, _ = build_glyphs(page, cls=c, size_range=band)
                for res in recognise(glyphs, runs, lib):
                    if not res['is_label']:
                        continue
                    P = [p for i in res['idx'] for pt in glyphs[i]['parts'] for p in pt['pts']]
                    x0, x1 = min(p[0] for p in P), max(p[0] for p in P)
                    y0, y1 = min(p[1] for p in P), max(p[1] for p in P)
                    gd = sum(glyphs[i]['d'] for i in res['idx']) / len(res['idx'])
                    # a decimal point touching the run (inside, or just beyond either end) => dimension
                    gh = max(y1 - y0, 1e-6)
                    if any(x0 - 1.2 * gd <= dx <= x1 + 1.2 * gd and y0 + 0.5 * gh <= dy <= y1 + 0.3 * gh
                           for dx, dy in dots):
                        continue
                    cx = (x0 + x1) / 2; cy = (y0 + y1) / 2
                    out.append({'text': res['text'], 'pt': (cx * M, (H - cy) * M),
                                'ambiguous': False})
                return out
            if cls == 'auto':
                bands = [b for b in size_groups(page) if b[1] >= 4.5]
                cand = {'size %.1f-%.1f' % b: _labels_for('upper', b) for b in bands}
                if len(bands) > 1:
                    merged = (min(b[0] for b in bands), max(b[1] for b in bands))
                    cand['size %.1f-%.1f (all text)' % merged] = _labels_for('upper', merged)
                if not cand:
                    cand = {'upper': _labels_for('upper')}
                def _score(ls):
                    # plot numbers come one per polygon; dimensions come several per polygon
                    cnt = collections.Counter()
                    for L in ls:
                        hit = [k for k, p in enumerate(polys) if _inpoly(p['rings'], L['pt'])]
                        if hit:
                            cnt[min(hit, key=lambda k: polys[k]['area'])] += 1
                    ones = sum(1 for v in cnt.values() if v == 1)
                    multi = sum(1 for v in cnt.values() if v > 1)
                    if ones == 0:
                        return (-10 ** 6, -len(ls))          # a band that labels nothing is never chosen
                    return (ones - multi, -len(ls))
                cls = max(cand, key=lambda c: _score(cand[c]))
                labels = cand[cls]
            else:
                labels = _labels_for(cls)

    # --- spatial join: label centroid inside smallest containing polygon ---
    assign = collections.defaultdict(list)
    unplaced = 0
    for li, L in enumerate(labels):
        hit = [k for k, p in enumerate(polys) if _inpoly(p['rings'], L['pt'])]
        if hit:
            assign[min(hit, key=lambda k: polys[k]['area'])].append(li)
        else:
            unplaced += 1

    def _poly_centroid(rings):
        r = rings[0][:-1]
        return (sum(p[0] for p in r) / len(r), sum(p[1] for p in r) / len(r))

    def _resolve(k, ls):
        """Pick ONE plot label for a polygon.

        On dense layout plans many dimension digits fall inside a plot, so
        concatenating them is wrong. Prefer the token nearest the polygon
        centroid; if several plausible tokens tie, mark the result for review
        rather than guessing. On sketches the glyph recogniser already returns
        whole labels, so a lone label passes through unchanged.
        """
        if not ls:
            return None, 0
        if len(ls) == 1:
            return labels[ls[0]]['text'], 1
        if kind == 'layout':
            # A named reserve (PARK/PP1/PP2/OSR/EWS) is unambiguous even amid
            # dimension noise; take it if exactly one is present. Otherwise the
            # numeric tokens here are plot dimensions, not the plot number, and
            # picking one would be a guess -- leave unlabelled (see the .md).
            named = [i for i in ls if not labels[i]['text'].replace('.', '').isdigit()]
            if len(named) == 1:
                return labels[named[0]]['text'], 1
            return None, len(ls)
        # sketch: several labels landed in one polygon (a missing dividing line, or a
        # dimension that slipped through). The plot number is the one nearest the
        # polygon centroid -- dimensions hug the edges; the rest go to alt_labels.
        c = _poly_centroid(polys[k]['rings'])
        order = sorted(ls, key=lambda i: math.dist(labels[i]['pt'], c))
        _resolve.alts[k] = '/'.join(labels[i]['text'] for i in order[1:])
        return labels[order[0]]['text'], len(ls)

    _resolve.alts = {}
    feats = []
    for k, p in enumerate(polys):
        ls = assign.get(k, [])
        txt, ncand = _resolve(k, ls)
        ring = p['rings'][0]
        per = sum(math.dist(ring[i], ring[i + 1]) for i in range(len(ring) - 1))
        feats.append({'type': 'Feature', 'id': k + 1, 'properties': {
            'poly_id': k + 1, 'survey_no': survey, 'plot_no': txt or None,
            'numbered': bool(txt),
            'area_sqm': round(p['area'], 3), 'area_are': round(p['area'] / 100, 4),
            'area_hect': round(p['area'] / 10000, 6), 'perimeter_m': round(per, 3),
            'n_vertices': len(ring) - 1, 'n_holes': len(p['rings']) - 1,
            'review_needed': bool(txt) and ncand > 1,
            'n_labels_in_polygon': ncand,
            'alt_labels': _resolve.alts.get(k) or None,
            'label_source': ('pdf_text_layer' if kind == 'layout'
                             else 'glyph_outline') if txt else 'none'},
            'geometry': {'type': 'Polygon', 'coordinates': p['rings']}})

    total = sum(f['properties']['area_sqm'] for f in feats)
    md = {'source_pdf': os.path.basename(pdf_path), 'document_kind': kind, 'text_class': cls,
          'village': village, 'taluk': taluk, 'survey_no': survey,
          'sheet_scale': f'1:{scale}', 'units': 'metres (ground)',
          'metres_per_pdf_point': round(M, 9),
          'snap_tolerance_m': round(tol * M, 4),
          'max_node_displacement_m': round(st['shift'] * M, 4),
          'arrows_dropped': n_arrows, 'slivers_merged': n_slivers,
          'stated_area_sqm': area, 'polygon_total_sqm': round(total, 2),
          'crs_note': ('Local planar CRS in ground metres, origin at page '
                       'bottom-left. NOT georeferenced -- apply a 2-point '
                       'Helmert fit to place on the ground.'),
          'geometry_note': ('Polygons are DERIVED: nodes snapped within the '
                            'tolerance above and dangling edges pruned. Use the '
                            'lines file for measurement-grade work.'),
          'labels_note': ('Plot numbers on dense LAYOUT plans are not extracted '
                          'automatically: dimension digits inside each plot are '
                          'indistinguishable from the plot number by position or '
                          'size, so only unambiguous named reserves (PARK/PP/OSR/'
                          'EWS) are labelled. Read the plot numbers off the sheet. '
                          'Geometry is unaffected. On FMB sketches with a glyph '
                          'library this note does not apply.') if kind=='layout' else
                         ('Plot numbers read from PDF vector glyph outlines via '
                          'shape matching; small rotated dimension text excluded.')}

    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(pdf_path))[0]
    crs = {'type': 'name', 'properties': {'name': 'urn:ogc:def:cs:local:sheet-metres'}}

    json.dump({'type': 'FeatureCollection', 'name': f'{stem}_parcels',
               'metadata': md, 'crs': crs, 'features': feats},
              open(f'{out_dir}/{stem}_parcels.geojson', 'w'), indent=1)

    lf = []
    for i, s in enumerate(segs):
        p0 = [round(s['p0'][0] * M, 4), round(s['p0'][1] * M, 4)]
        p1 = [round(s['p1'][0] * M, 4), round(s['p1'][1] * M, 4)]
        lf.append({'type': 'Feature', 'id': i + 1, 'properties': {
            'id': i + 1, 'survey_no': survey, 'layer': s['layer'],
            'linewidth_pt': s['lw'], 'dashed': s['dashed'],
            'length_m': round(math.dist(p0, p1), 4)},
            'geometry': {'type': 'LineString', 'coordinates': [p0, p1]}})
    json.dump({'type': 'FeatureCollection', 'name': f'{stem}_lines',
               'metadata': md, 'crs': crs, 'features': lf},
              open(f'{out_dir}/{stem}_lines.geojson', 'w'), indent=1)

    _fields = ['poly_id', 'survey_no', 'plot_no', 'numbered', 'area_sqm', 'area_are', 'area_hect',
               'perimeter_m', 'n_vertices', 'n_holes', 'review_needed', 'n_labels_in_polygon', 'alt_labels', 'label_source']
    with open(f'{out_dir}/{stem}_parcels.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=_fields)
        w.writeheader()
        for f in sorted(feats, key=lambda f: (f['properties']['plot_no'] or 'zzz')):
            w.writerow(f['properties'])

    if feats:
        _preview(feats, md, f'{out_dir}/{stem}_preview.png')

    numbered = sum(1 for f in feats if f['properties']['plot_no'])
    dev = (f'{100 * (total - area) / area:+.2f}%' if area else 'n/a')
    print(f'{stem}: {kind} | polygons {len(feats)} | numbered {numbered} | '
          f'unplaced labels {unplaced}')
    print(f'  total {total:.1f} m2 | stated {area} | deviation {dev}')
    print(f'  scale 1:{scale} ({M:.6f} m/pt) | max node shift '
          f"{st['shift'] * M:.3f} m (tol {tol * M:.2f} m) | components {st['components']}")
    print(f'  wrote {out_dir}/{stem}_parcels.geojson (+ .csv, _lines.geojson, _preview.png)')
    run.last_stats = {'polygons': len(feats), 'numbered': numbered, 'unplaced_labels': unplaced,
                      'arrows_dropped': n_arrows, 'slivers_merged': n_slivers,
                      'total_sqm': round(total, 2), 'max_node_shift_m': round(st['shift'] * M, 4),
                      'components': st['components'], 'text_class': cls,
                      'review': sum(1 for f in feats if f['properties']['review_needed'])}
    return feats


def _preview(feats, md, path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.path import Path
    except Exception:
        return
    n = len(feats)
    fig, ax = plt.subplots(figsize=(max(11, min(26, 9 + n / 14)),) * 2)
    for f in feats:
        rs = f['geometry']['coordinates']; v = []; c = []
        for r in rs:
            v += r; c += [Path.MOVETO] + [Path.LINETO] * (len(r) - 2) + [Path.CLOSEPOLY]
        lab = f['properties']['plot_no']
        ax.add_patch(mpatches.PathPatch(Path(v, c),
                     fc='#ffe9b0' if lab else '#c9c9c9', ec='k', lw=.5))
        r0 = rs[0][:-1]
        cx = sum(p[0] for p in r0) / len(r0); cy = sum(p[1] for p in r0) / len(r0)
        ax.text(cx, cy, lab or '', ha='center', va='center',
                fontsize=max(2.6, min(7, 260 / max(12, n))), color='b')
    ax.autoscale_view(); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f"{md.get('survey_no') or md['source_pdf']} - {n} polygons, "
                 f"{sum(1 for f in feats if f['properties']['plot_no'])} numbered "
                 f"({md['sheet_scale']})")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def show_titleblock(pdf_path):
    """Print the words in the top strip so you can read scale/area/village."""
    import pdfplumber
    page = pdfplumber.open(pdf_path).pages[0]
    H = float(page.height)
    top = [w for w in page.extract_words() if H - w['top'] > H - H * 0.13]
    for w in sorted(top, key=lambda w: (round(w['top']), w['x0'])):
        print(f"  {w['text']}")


def main():
    ap = argparse.ArgumentParser(description='FMB / CMDA cadastral PDF -> GeoJSON')
    ap.add_argument('pdf')
    ap.add_argument('--scale', type=float, help='sheet scale denominator, e.g. 500 for 1:500')
    ap.add_argument('--area', type=float, help='stated area in m2 (for a closure check)')
    ap.add_argument('--kind', choices=['auto', 'sketch', 'layout'], default='auto')
    ap.add_argument('--glyphlib', help='glyph exemplar JSON (FMB sketches only)')
    ap.add_argument('--village'); ap.add_argument('--taluk'); ap.add_argument('--survey')
    ap.add_argument('--clip', nargs=4, type=float, metavar=('X0', 'X1', 'Y0', 'Y1'),
                    help='keep only geometry inside this box (PDF pt, y-up)')
    ap.add_argument('--snap', type=float, default=2.7,
                    help='node-snap tolerance in PDF points at A0 (scaled by sheet size)')
    ap.add_argument('--cls', choices=['auto', 'upper', 'lower'], default='auto',
                    help='which text size class holds the plot numbers')
    ap.add_argument('-o', '--out', default='out')
    ap.add_argument('--show-titleblock', action='store_true',
                    help='print the title-block words and exit')
    a = ap.parse_args()

    if a.show_titleblock:
        show_titleblock(a.pdf); return
    if not a.scale:
        ap.error('--scale is required (read it off the title block, or use --show-titleblock)')
    run(a.pdf, a.scale, a.out, kind=a.kind, glyphlib=a.glyphlib, area=a.area,
        village=a.village, taluk=a.taluk, survey=a.survey,
        clip=tuple(a.clip) if a.clip else None, snap_pt=a.snap, cls=a.cls)


if __name__ == '__main__':
    main()
