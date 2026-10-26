"""Tidy a hand-georeferenced village so it reads neat and clean, WITHOUT moving any parcel.

    python -m autogeoref.village_clean 35_04_074 Thirukatchur --src <geojson> --out <geojson>

Akash, 2026-09-24: internal lines were smudged (double lines, zigzags, needles, slivers). A first
attempt rebuilt every survey from its PDF at a best-fit pose; that moved parcels off the place
he had put them (47B plot 1 by 7.7 m, 48A plot 2 shrunk to two thirds) and he rejected it: the
placement is his and is right. So his placement is kept and only the drawing is cleaned:

  1. OUTLINE: each survey keeps its own outline; needles and notches narrower than 2 x TIDY_R
     are smoothed out of it, piece by piece, and only when the piece is small.
  2. INTERIOR: the survey's plots are made to tile that outline exactly (no overlap, no gap,
     one shared line between two plots). The PDF's own straight internal lines are used ONLY
     when every plot stays within PLOT_IOU_MIN / PLOT_AREA_TOL of his plot of the same number;
     otherwise his own plots are kept and tidied.
  3. SEAMS, railway land first then outward: where a neighbour's line is within SNAP_M the two
     lines become one (outer points move onto the neighbour, never internal points). Empty
     ground a thin strip wide between parcels is given to the side it changes least.
  4. The report measures every survey and plot against the placement it was given (overlap, area, largest
     boundary move), not only against the PDF, so a reduced or displaced parcel cannot pass.
"""
import argparse
import logging
import re
import sys
import zipfile

import geopandas as gpd
import numpy as np
from scipy.optimize import minimize
from shapely import affinity, set_precision
from shapely.geometry import LineString, MultiPoint, Point, Polygon
from shapely.ops import nearest_points, substring, unary_union
from shapely.validation import make_valid

from . import engine, fmb_finish, fmb_on_base, neighbours, paths, shift_puvi

log = logging.getLogger(__name__)

UTM = 32644
ACRE = 4046.8564224
SNAP_TOL_M = 6.0        # a seam wider than this is left open and reported, never forced
VERTEX_SNAP_M = 1.0     # within this of a neighbour's vertex, share that vertex
FRAGMENT_M2 = 0.5
GRID = 0.001            # final precision grid (1 mm): shared vertices stay shared


# --------------------------------------------------------------------------- geometry helpers

def _polys(g):
    return shift_puvi._polygons_of(g) if g is not None and not g.is_empty else []


def _clean_poly(g):
    """Valid polygonal geometry, fragments dropped; None when nothing is left."""
    if g is None or g.is_empty:
        return None
    if not g.is_valid:
        g = make_valid(g)
    parts = [p for p in _polys(g) if p.area >= FRAGMENT_M2]
    if not parts:
        return None
    out = unary_union(parts)
    return out if out.is_valid else out.buffer(0)


def _iou(a, b):
    i = a.intersection(b).area
    u = a.area + b.area - i
    return i / u if u else 0.0


def _key(c):
    return (round(c[0] / 0.01), round(c[1] / 0.01))


def _rebuild(geom, table):
    """Rebuild a polygon replacing any vertex found in `table` (1 cm keys)."""
    def ring(coords):
        return [table.get(_key(c), (float(c[0]), float(c[1]))) for c in coords]
    parts = []
    for p in _polys(geom):
        parts.append(Polygon(ring(p.exterior.coords), [ring(r.coords) for r in p.interiors]))
    return _clean_poly(unary_union(parts) if len(parts) > 1 else (parts[0] if parts else None))


# --------------------------------------------------------------------------- 1. regenerate

def _apply(geom, pose):
    ang, tx, ty, ox, oy = pose
    return affinity.rotate(affinity.translate(geom, tx, ty), ang, origin=(ox, oy))


def best_pose(sheet_gdf, current):
    """Rigid pose putting the sheet over the survey's current plots. current: [(label, geom)]."""
    outline = unary_union([fmb_on_base._valid(g) for g in sheet_gdf.geometry])
    target = unary_union([g for _l, g in current])
    ang0, shift0, _s = shift_puvi.rigid_fit(outline, target)
    ox, oy = target.centroid.x, target.centroid.y
    cur_by_label = {}
    for lab, g in current:
        if lab:
            cur_by_label.setdefault(lab, []).append(g)
    cur_by_label = {k: unary_union(v) for k, v in cur_by_label.items()}
    labelled = [(str(l), fmb_on_base._valid(g)) for l, g in zip(sheet_gdf["plot_no"], sheet_gdf.geometry)
                if str(l) in cur_by_label]
    lab_area = sum(g.area for _l, g in labelled) or 1.0

    def score(pose):
        s = _iou(_apply(outline, pose), target)
        if labelled:
            s += sum(_apply(g, pose).intersection(cur_by_label[l]).area for l, g in labelled) / lab_area
        return s

    best = None
    for flip in (0.0, 180.0):
        start = (ang0 + flip, float(shift0[0]), float(shift0[1]))
        f = lambda x: -score((x[0], x[1], x[2], ox, oy))
        res = minimize(f, np.array(start), method="Nelder-Mead",
                       options={"initial_simplex": np.array([start,
                                                              [start[0] + 1.0, start[1], start[2]],
                                                              [start[0], start[1] + 1.0, start[2]],
                                                              [start[0], start[1], start[2] + 1.0]]),
                                "xatol": 1e-3, "fatol": 1e-6, "maxiter": 600})
        if best is None or -res.fun > best[0]:
            best = (-res.fun, (res.x[0], res.x[1], res.x[2], ox, oy))
    return best[1]


def regenerate(village, survey, current):
    """[(plot_no, geom)] from the PDF sheet at the pose fitted to `current`; None without a sheet."""
    sheet_gdf, outline = fmb_on_base.sheet_of(village, survey)
    if outline is None:
        return None, None
    sheet_gdf = sheet_gdf.copy()
    sheet_gdf["plot_no"] = ["" if str(x).strip().lower() in ("", "nan", "none", "0") else str(x).strip()
                            for x in sheet_gdf.get("plot_no", [""] * len(sheet_gdf))]
    pose = best_pose(sheet_gdf, current)
    plots = [(lab, _clean_poly(_apply(fmb_on_base._valid(g), pose)))
             for lab, g in zip(sheet_gdf["plot_no"], sheet_gdf.geometry)]
    return [(l, g) for l, g in plots if g is not None], _apply(outline, pose)


# --------------------------------------------------------------------------- 1b. tidy in place

TIDY_R = 0.75           # needles and notches narrower than 2 x this are smoothed out of an outline
TIDY_PIECE_M2 = 12.0    # ...one piece at a time, and only a piece this small (a real bay stays)
SNAP_M = 1.5            # a neighbour's line this close is the same line: the two become one
PLOT_IOU_MIN = 0.90     # the PDF's internal lines are used only when every plot stays this close
PLOT_AREA_TOL = 0.10    # ...to his plot of the same number, and within this share of its area
_MITRE = dict(join_style=2, mitre_limit=50.0)


def _small(piece):
    w = 2.0 * piece.area / piece.length if piece.length else 0.0
    return piece.area <= TIDY_PIECE_M2 and w <= 2.0 * TIDY_R


def tidy_outline(body):
    """The survey's outline without needles, notches and zigzags narrower than 2 x TIDY_R.

    Mitred opening and closing put every straight edge and every corner back exactly where it
    was (a high mitre limit, so even a strip's 10-degree tip survives); what does not come back
    is a feature too thin to have been drawn on purpose. Each removed or added piece must be
    small, so a narrow lane or a real bay in the outline is never touched.
    """
    body = _clean_poly(body)
    if body is None:
        return None
    opened = body.buffer(-TIDY_R, **_MITRE).buffer(TIDY_R, **_MITRE)
    cut = [q for q in _polys(body.difference(opened)) if q.area > 1e-4 and _small(q)]
    out = _clean_poly(body.difference(unary_union(cut))) if cut else body
    out = out if out is not None else body
    closed = out.buffer(TIDY_R, **_MITRE).buffer(-TIDY_R, **_MITRE)
    add = [q for q in _polys(closed.difference(out)) if q.area > 1e-4 and _small(q)]
    if add:
        grown = _clean_poly(unary_union([out] + add))
        out = grown if grown is not None else out
    return out


def partition(plots, outline):
    """Plots that tile `outline` exactly. Labels never change and no plot is ever emptied.

    Every plot is clipped to the outline; where two plots overlap the smaller one yields the
    hairline; a detached crumb of a plot is handed back; ground no plot covers goes to the plot
    bordering it, split along the straight extensions of the internal lines so that each line
    runs on to the outline instead of stepping.
    """
    out = []
    for lab, g in plots:
        g2 = _clean_poly(g.intersection(outline)) if g is not None else None
        out.append([lab, g2])
    taken = None
    for i in sorted(range(len(out)), key=lambda i: -(out[i][1].area if out[i][1] is not None else 0.0)):
        g = out[i][1]
        if g is None:
            continue
        if taken is not None and g.intersects(taken):
            g2 = _clean_poly(g.difference(taken))
            if g2 is not None:
                g = g2
        parts = _polys(g)
        if len(parts) > 1:                       # keep the plot's own body, hand crumbs back
            big = max(parts, key=lambda q: q.area)
            g = unary_union([q for q in parts if q is big or q.area >= 0.25 * big.area])
        out[i][1] = g
        taken = g if taken is None else unary_union([taken, g])
    plots2 = [(l, g) for l, g in out if g is not None]
    body = unary_union([g for _l, g in plots2])
    rest = outline.difference(body)
    for piece in sorted([q for q in _polys(rest) if q.area > 1e-4], key=lambda q: -q.area):
        body = unary_union([g for _l, g in plots2])
        plots2 = _give(plots2, piece, body)
    return [(l, _clean_poly(g) or g) for l, g in plots2]


def _outer(g):
    """The outer lines of a (multi)polygon, holes left out."""
    from shapely.geometry import MultiLineString
    return MultiLineString([list(p.exterior.coords) for p in _polys(g)])


def _match(new, cur):
    """For each plot of `new`, (label, IoU, area ratio) against his best plot of the same number."""
    by = {}
    for lab, g in cur:
        by.setdefault(lab, []).append(g)
    res = []
    for lab, g in new:
        cands = by.get(lab, [])
        if not cands:
            res.append((lab, 0.0, 0.0))
            continue
        c = max(cands, key=lambda c: _iou(g, c))
        res.append((lab, _iou(g, c), g.area / c.area if c.area else 0.0))
    return res


def interior_gate(new, cur):
    """True when the PDF's interior changes none of his plots visibly."""
    if len(new) != len(cur) or sorted(l for l, _g in new) != sorted(l for l, _g in cur):
        return False
    return all(i >= PLOT_IOU_MIN and abs(r - 1.0) <= PLOT_AREA_TOL for _l, i, r in _match(new, cur))


# --------------------------------------------------------------------------- 2. conflate

PARALLEL_DEG = 25.0     # an edge "runs along" a neighbour when their directions agree this well
WIDTH_SHARE = 0.35      # no point moves more than this share of the survey's own mean width


def _bearing(a, b):
    return np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0])) % 180.0


def _tangent(rings, q):
    """Undirected bearing of the neighbour boundary at q."""
    best = None
    for R in rings:
        d = R.distance(q)
        if best is None or d < best[0]:
            best = (d, R)
    R = best[1]
    s = R.project(q)
    a = R.interpolate(max(0.0, s - 0.5))
    b = R.interpolate(min(R.length, s + 0.5))
    return _bearing((a.x, a.y), (b.x, b.y))


def _angle_diff(x, y):
    d = abs(x - y) % 180.0
    return min(d, 180.0 - d)


def snap_outer_nodes(plots, seniors, tol=SNAP_TOL_M):
    """Move OUTER points that sit on a real seam onto the finished neighbours' boundary.

    A point moves only when one of its outer edges runs along the neighbour's edge (both ends
    within reach, directions within PARALLEL_DEG) - a corner that merely comes near a neighbour
    stays. And it never moves more than WIDTH_SHARE of the survey's mean width, so a small or
    narrow survey cannot be collapsed onto its neighbours (41 lost everything on 2026-09-24).
    """
    body = unary_union([g for _l, g in plots])
    sb = seniors.boundary
    rings = [LineString(r.coords) for p in _polys(seniors) for r in [p.exterior, *p.interiors]]
    svert = MultiPoint([c for R in rings for c in R.coords])
    mean_w = 2.0 * body.area / body.length if body.length else 0.0
    cap = min(tol, WIDTH_SHARE * mean_w) if mean_w else tol
    nodes = shift_puvi._nodes([g for _l, g in plots], grid=0.01)
    # outer ring order, to know each outer point's two outer edges
    nbrs = {}
    for poly in _polys(body):
        cs = list(poly.exterior.coords)[:-1]
        n = len(cs)
        for i in range(n):
            k = _key(cs[i])
            nbrs.setdefault(k, []).extend([cs[i - 1], cs[(i + 1) % n]])
    table, moves = {}, []
    for k, c in nodes.items():
        if k not in nbrs:
            continue                                   # an internal point: never moves
        pt = Point(c)
        d = sb.distance(pt)
        if d > cap or d < 1e-4:
            continue
        inside = seniors.contains(pt)
        q = nearest_points(sb, pt)[0]
        runs_along = False
        for other in nbrs[k]:
            if sb.distance(Point(other)) > tol:
                continue
            # the neighbour's direction is read at the middle of the edge, never at a corner
            mid = Point((c[0] + other[0]) / 2.0, (c[1] + other[1]) / 2.0)
            if sb.distance(mid) > tol:
                continue
            tan = _tangent(rings, nearest_points(sb, mid)[0])
            if _angle_diff(_bearing(c, other), tan) <= PARALLEL_DEG:
                runs_along = True
                break
        if not runs_along:
            continue
        v = nearest_points(svert, pt)[0]
        if pt.distance(v) <= min(VERTEX_SNAP_M, cap) and pt.distance(v) <= d + 0.5:
            q = v
        if not inside:
            seg = LineString([c, (q.x, q.y)])
            if seg.intersection(body).length > 0.25 * seg.length + 0.05:
                continue                               # the path runs back across the survey
        table[k] = (q.x, q.y)
        moves.append(pt.distance(q))
    out = []
    for lab, g in plots:
        ng = _rebuild(g, table)
        out.append((lab, ng if ng is not None else g))
    return out, moves


def follow_seams(plots, seniors, tol=SNAP_TOL_M):
    """Replace every edge lying on a finished neighbour's boundary with that boundary itself."""
    rings = [LineString(r.coords) for p in _polys(seniors) for r in [p.exterior, *p.interiors]]

    def ring_of(c):
        for i, R in enumerate(rings):
            if R.distance(Point(c)) < 1e-3:
                return i
        return None

    def path(R, a, b):
        L, da, db = R.length, R.project(Point(a)), R.project(Point(b))
        if abs(da - db) < 1e-6:
            return None
        if da < db:
            fwd = list(substring(R, da, db).coords)
            back = list(substring(R, db, L).coords) + list(substring(R, 0, da).coords)[1:]
            back = back[::-1]
        else:
            fwd = list(substring(R, da, L).coords) + list(substring(R, 0, db).coords)[1:]
            back = list(substring(R, db, da).coords)[::-1]
        seg = LineString([a, b])
        best = None
        for cand in (fwd, back):
            if len(cand) < 2:
                continue
            ln = LineString(cand)
            h = ln.hausdorff_distance(seg)
            if h <= tol and ln.length <= 1.5 * seg.length + 1.0 and (best is None or h < best[0]):
                best = (h, cand)
        return best[1] if best else None

    def ring_new(coords):
        cs = list(coords)[:-1]
        n = len(cs)
        out = []
        for i in range(n):
            a, b = cs[i], cs[(i + 1) % n]
            out.append(a)
            ra, rb = ring_of(a), ring_of(b)
            if ra is not None and ra == rb:
                mid = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)
                if rings[ra].distance(Point(mid)) <= tol:
                    p = path(rings[ra], a, b)
                    if p:
                        out.extend(p[1:-1])
        return out

    res = []
    for lab, g in plots:
        parts = []
        for p in _polys(g):
            parts.append(Polygon(ring_new(p.exterior.coords), [r.coords for r in p.interiors]))
        ng = _clean_poly(unary_union(parts) if parts else None)
        res.append((lab, ng if ng is not None else g))
    return res


MAX_LOSS = 0.10         # a conflation that costs a survey more than this of its drawing is retried gentler
MAX_PLOT_LOSS = 0.40


def _conflate_once(plots, seniors, tol):
    plots2, moves = snap_outer_nodes(plots, seniors, tol)
    plots2 = follow_seams(plots2, seniors, min(tol, 2.0))
    out = []
    for lab, g in plots2:
        ng = _clean_poly(g.difference(seniors)) if g.intersects(seniors) else g
        out.append((lab, ng))
    return out, moves


def conflate(plots, seniors, tol=SNAP_TOL_M):
    """Snap, follow, cut - retried with a smaller reach whenever it would damage the survey."""
    if seniors is None or seniors.is_empty:
        return plots, []
    before = {i: g.area for i, (_l, g) in enumerate(plots)}
    total = sum(before.values())
    for t_ in (tol, tol / 2.0, 1.5, 0.0):
        if t_ > 0:
            cand, moves = _conflate_once(plots, seniors, t_)
        else:                                         # last resort: only cut what overlaps
            cand = [(l, _clean_poly(g.difference(seniors)) if g.intersects(seniors) else g) for l, g in plots]
            moves = []
        lost_total = 1.0 - sum(g.area for _l, g in cand if g is not None) / total if total else 0.0
        worst_plot = max((1.0 - (cand[i][1].area if cand[i][1] is not None else 0.0) / a)
                         for i, a in before.items() if a > 0)
        if (lost_total <= MAX_LOSS and worst_plot <= MAX_PLOT_LOSS) or t_ == 0.0:
            return [(l, g) for l, g in cand if g is not None], moves
    return plots, []


# --------------------------------------------------------------------------- 3. close seams

SEAM_MAX_M = 6.0        # a gap between PDF neighbours up to this is a seam and is closed
SLIVER_MAX_M = 1.5      # between ANY two parcels a gap this thin is never a road: closed too
HOLE_MAX_M2 = 300.0     # an enclosed hole this small is a seam, not unplaced ground
SEAM_ROUNDS = 4
SEAM_SHARE_MAX = 0.15   # one side takes a whole seam only if its plot grows by no more than this


def _internal_edges(plots, body):
    """Edges of the survey's plots that are NOT on its outline: its internal lines."""
    out = []
    ext = body.boundary
    for _l, g in plots:
        for p in _polys(g):
            cs = list(p.exterior.coords)
            for i in range(len(cs) - 1):
                a, b = cs[i], cs[i + 1]
                if ext.distance(Point((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)) > 0.01:
                    out.append((a, b))
    return out


def _extension_lines(plots, body, region, reach):
    """Straight continuations of internal lines that end on the outline next to `region`."""
    lines = []
    for a, b in _internal_edges(plots, body):
        for inner, end in ((a, b), (b, a)):
            pe = Point(end)
            if body.boundary.distance(pe) > 0.01 or region.distance(pe) > 0.05:
                continue
            d = np.array(end[:2]) - np.array(inner[:2])
            n = float(np.hypot(*d))
            if n < 1e-6:
                continue
            d = d / n
            far = (end[0] + d[0] * (reach + 3.0), end[1] + d[1] * (reach + 3.0))
            lines.append(LineString([end[:2], far]))
    return lines


def _gap_between(ma, mb, width):
    """The ground between two parcels' facing edges (morphological closing of their union)."""
    r = width / 2.0 + 0.6
    u = unary_union([ma, mb])
    closed = u.buffer(r, join_style=2, mitre_limit=2.0).buffer(-r, join_style=2, mitre_limit=2.0)
    g = closed.difference(u)
    return [q for q in _polys(g) if q.area > 0.05 and q.distance(ma) < 0.02 and q.distance(mb) < 0.02]


def _give(plots, piece, body):
    """Merge `piece` into the survey's plots, each part to the plot it lies against.

    The piece is first cut along the straight continuations of the internal lines that end on
    it, reaching only across the piece (a whole-survey reach once cut a thin strip in the wrong
    places and handed 42B's corner wedge to a plot 128 m away). A part that still lies against
    several plots is shared out to the nearest of them, never given whole to one.
    """
    from shapely.ops import split as _split
    pieces = [piece]
    mean_w = 2.0 * piece.area / piece.length if piece.length else 0.0
    span = max(piece.bounds[2] - piece.bounds[0], piece.bounds[3] - piece.bounds[1])
    reach = min(span, 4.0 * mean_w + 2.0)
    for ln in _extension_lines(plots, body, piece, reach):
        nxt = []
        for q in pieces:
            try:
                parts = [x for x in _split(q, ln).geoms if x.area > 1e-4]
            except Exception:
                parts = [q]
            nxt.extend(parts if parts else [q])
        pieces = nxt
    plots = list(plots)
    for q in pieces:
        qb = q.buffer(0.02)
        contact = {i: g.boundary.intersection(qb).length for i, (_l, g) in enumerate(plots)}
        contact = {i: L for i, L in contact.items() if L > 0.05}
        if not contact:
            continue
        top = max(contact, key=contact.get)
        if len(contact) == 1 or contact[top] >= 0.8 * sum(contact.values()):
            gifts = {top: q}
        else:
            gifts = _split_nearest(q, {i: plots[i][1] for i in contact})
        for i, part in gifts.items():
            lab, g = plots[i]
            merged = unary_union([g, part])
            merged = merged if merged.is_valid else merged.buffer(0)
            plots[i] = (lab, merged)
    return plots


def _split_nearest(q, owners):
    """{name: piece of q nearest that owner}: the gap cut down its middle line.

    Built from a Voronoi diagram of points 0.25 m apart along each owner's boundary next to the
    gap, so each side moves at most half the gap and a long sliver along several parcels is
    shared out along its length instead of growing one parcel a tail past the others.
    """
    zone = q.buffer(0.5)
    pts, who, seen = [], [], set()
    for name, b in owners.items():
        edge = b.boundary.intersection(zone)
        for ln in (getattr(edge, "geoms", None) or [edge]):
            if ln.is_empty or ln.geom_type not in ("LineString", "LinearRing") or ln.length == 0:
                continue
            n = max(2, int(ln.length / 0.25) + 1)
            for i in range(n):
                c = ln.interpolate(i * ln.length / (n - 1))
                k = (round(c.x, 3), round(c.y, 3))
                if k in seen or any(o != name and ob.boundary.distance(c) < 0.005 for o, ob in owners.items()):
                    continue                     # a point two owners share decides nothing
                seen.add(k)
                pts.append(c)
                who.append(name)
    if len(set(who)) < 2:
        return {who[0]: q} if who else {}
    from shapely import voronoi_polygons
    cells = voronoi_polygons(MultiPoint(pts), extend_to=q.envelope.buffer(5.0), ordered=True)
    by = {}
    for name, cell in zip(who, cells.geoms):
        by.setdefault(name, []).append(cell)
    out = {}
    for name, cs in by.items():
        piece = unary_union(cs).intersection(q)
        piece = unary_union([x for x in _polys(piece) if x.area > 0.01]) if not piece.is_empty else piece
        if not piece.is_empty:
            out[name] = piece
    return out


def _receiving_plot_area(plots, piece):
    best, blen = None, 0.0
    qb = piece.buffer(0.02)
    for _l, g in plots:
        L = g.boundary.intersection(qb).length
        if L > blen:
            best, blen = g, L
    return best.area if best is not None else 0.0


def close_seams(village, done):
    """Fill the EMPTY ground between parcels that the eye reads as a gap. Returns notes.

    Pairwise distance was the wrong yardstick: two PDF neighbours can be apart because a third
    parcel lies between them, which looks perfectly clean. What looks wrong is empty ground:
    a thin strip or small hole with parcels on both sides. Found by morphological closing of the
    fabric; each piece goes to the bordering survey for which it is the smallest change, split
    along that survey's internal lines. Between PDF neighbours a gap up to SEAM_MAX_M wide is
    closed; between parcels the sheets do NOT print as neighbours only up to SLIVER_MAX_M, so a
    real lane between them is never filled in.
    """
    notes = []
    r = SEAM_MAX_M / 2.0
    for _round in range(SEAM_ROUNDS):
        bodies = {s: unary_union([g for _l, g in p]) for s, p in done.items() if p}
        adj = set()
        for s in bodies:
            for n in neighbours.printed(village, s):
                if n in bodies and n != s:
                    adj.add(tuple(sorted((s, n))))
        u = unary_union(list(bodies.values()))
        closed = u.buffer(r, join_style=2, mitre_limit=2.0).buffer(-r, join_style=2, mitre_limit=2.0)
        pieces = [q for q in _polys(closed.difference(u)) if q.area > 0.01]
        for poly in _polys(u):
            for ring in poly.interiors:
                h = Polygon(ring)
                if 0.01 < h.area <= HOLE_MAX_M2 and not any(h.equals(q) for q in pieces):
                    pieces.append(h)
        changed = 0
        for q in sorted(pieces, key=lambda x: x.area):
            bodies = {s: unary_union([g for _l, g in p]) for s, p in done.items() if p}
            around = [s for s in bodies if bodies[s].distance(q) < 0.02 and
                      bodies[s].boundary.intersection(q.buffer(0.02)).length > 0.05]
            if len(around) < 2:
                continue                                  # a notch in one parcel's own outline
            q2 = q.difference(unary_union([bodies[s] for s in around]))
            q2 = unary_union([x for x in _polys(q2) if x.area > 0.01]) if not q2.is_empty else q2
            if q2.is_empty:
                continue
            width = 2.0 * q2.area / q2.length if q2.length else 0.0
            is_seam = any(tuple(sorted((a, b))) in adj for a in around for b in around if a < b)
            if width > (SEAM_MAX_M if is_seam else SLIVER_MAX_M):
                continue
            lens = {s: bodies[s].boundary.intersection(q2.buffer(0.02)).length for s in around}
            # the SIDES of the gap; a parcel that only closes off its end gets none of it (a
            # railway strip touching the end of a seam grew a triangular tooth on 2026-09-24)
            sides = [s for s in around if lens[s] >= 0.2 * max(lens.values())]
            if len(sides) == 1 and width > SLIVER_MAX_M:
                continue
            one_seam = len(sides) <= 2 and min(lens[s] for s in sides) >= 0.5 * max(lens.values())
            share = {s: q2.area / max(_receiving_plot_area(done[s], q2), 1e-6) for s in sides}
            tgt = min(sides, key=lambda s: share[s])
            if one_seam and share[tgt] <= SEAM_SHARE_MAX:
                # one seam between two parcels: the side for which it is the smaller change takes
                # it, judged by the PLOT that receives it, so a narrow strip keeps its drawn width;
                # it takes only the part lying between the two, never a tail past its own end
                reach = 2.5 * width + 0.05
                between = q2.intersection(bodies[tgt].buffer(reach, join_style=2, mitre_limit=2.0))
                between = unary_union([x for x in _polys(between) if x.area > 0.01]) if not between.is_empty else between
                gifts = {tgt: between} if not between.is_empty else {}
                beyond = q2.difference(between) if not between.is_empty else q2
                for x in _polys(beyond):
                    if x.area > 0.01:
                        for s2, piece in _split_nearest(x, {s: bodies[s] for s in sides}).items():
                            gifts[s2] = unary_union([gifts[s2], piece]) if s2 in gifts else piece
            else:
                # several parcels along it, or too much for either side alone: down the middle
                gifts = _split_nearest(q2, {s: bodies[s] for s in sides})
            for tgt, piece in gifts.items():
                done[tgt] = _give(done[tgt], piece, bodies[tgt])
            changed += 1
            notes.append("%.1f m2 (%.2f m wide) between %s into %s" % (
                q2.area, width, "/".join(sorted(around)),
                "+".join("%s %.1f" % (s, gifts[s].area) for s in sorted(gifts))))
        if not changed:
            break
    return notes


# --------------------------------------------------------------------------- 4. finish the coverage

MERGE_M = 0.05          # two vertices this close are the same drawn point
SLIVER_PLOT_M = 0.5     # an UNNUMBERED plot thinner than this is a double line, not a plot


def absorb_slivers(plots):
    """Merge unnumbered hairline plots (a line drawn twice on the sheet) into their neighbour."""
    keep = [(l, g) for l, g in plots]
    changed = True
    while changed:
        changed = False
        for i, (l, g) in enumerate(keep):
            if l or g.length == 0 or 2.0 * g.area / g.length >= SLIVER_PLOT_M or len(keep) < 2:
                continue
            gb = g.buffer(0.02)
            j = max((k for k in range(len(keep)) if k != i),
                    key=lambda k: keep[k][1].boundary.intersection(gb).length)
            if keep[j][1].boundary.intersection(gb).length <= 0:
                continue
            keep[j] = (keep[j][0], unary_union([keep[j][1], g]))
            del keep[i]
            changed = True
            break
    return keep


def finish_coverage(done):
    """One coordinate per drawn point across the whole village, and nothing left in between.

    Vertices closer than MERGE_M become one (a zigzag of a few centimetres is gone and two
    lines meeting at a point meet at the SAME point); a vertex lying on a neighbour's edge is
    added to that edge (so a seam is one line in any GIS checker, not two that merely touch);
    pin-holes and hairline overlaps left by the earlier steps are given to one side.
    """
    from scipy.spatial import cKDTree
    items = [(s, i) for s in done for i in range(len(done[s]))]
    coords = []
    for s, i in items:
        for p in _polys(done[s][i][1]):
            for r in [p.exterior, *p.interiors]:
                coords.extend((float(x), float(y)) for x, y in list(r.coords)[:-1])
    if not coords:
        return done
    arr = np.array(coords)
    parent = list(range(len(arr)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in cKDTree(arr).query_pairs(MERGE_M):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    groups = {}
    for k in range(len(arr)):
        groups.setdefault(find(k), []).append(k)
    table = {}
    for members in groups.values():
        if len(members) < 2:
            continue
        cnt = {}
        for k in members:
            c = (arr[k][0], arr[k][1])
            cnt[c] = cnt.get(c, 0) + 1
        rep = max(cnt, key=lambda c: cnt[c])
        for k in members:
            table[(arr[k][0], arr[k][1])] = rep

    def ring(cs):
        out = []
        for c in list(cs)[:-1]:
            q = table.get((float(c[0]), float(c[1])), (float(c[0]), float(c[1])))
            if not out or (abs(q[0] - out[-1][0]) > 1e-9 or abs(q[1] - out[-1][1]) > 1e-9):
                out.append(q)
        while len(out) > 1 and out[0] == out[-1]:
            out.pop()
        return out if len(out) >= 3 else None

    for s, i in items:
        lab, g = done[s][i]
        parts = []
        for p in _polys(g):
            ext = ring(p.exterior.coords)
            if ext is None:
                continue
            holes = [h for h in (ring(r.coords) for r in p.interiors) if h is not None]
            parts.append(Polygon(ext, holes))
        ng = _clean_poly(unary_union(parts) if len(parts) > 1 else (parts[0] if parts else None))
        done[s][i] = (lab, ng if ng is not None else g)
    _onto_edges(done)
    return _faces(done)


def _onto_edges(done):
    """A vertex within MERGE_M of another plot's edge moves onto that edge (every copy of it)."""
    import shapely
    items = [(s, i) for s in done for i in range(len(done[s]))]
    segs, owner = [], []
    for k, (s, i) in enumerate(items):
        for p in _polys(done[s][i][1]):
            for r in [p.exterior, *p.interiors]:
                cs = list(r.coords)
                for a, b in zip(cs[:-1], cs[1:]):
                    segs.append(LineString([a, b]))
                    owner.append(k)
    if not segs:
        return done
    tree = shapely.STRtree(segs)
    table = {}
    for k, (s, i) in enumerate(items):
        for p in _polys(done[s][i][1]):
            for r in [p.exterior, *p.interiors]:
                for c in list(r.coords)[:-1]:
                    c = (float(c[0]), float(c[1]))
                    if c in table:
                        continue
                    v = Point(c)
                    best = None
                    for j in tree.query(v.buffer(MERGE_M)):
                        if owner[j] == k:
                            continue
                        sg = segs[j]
                        d = sg.distance(v)
                        if d < 1e-9 or d >= MERGE_M:
                            continue
                        e0, e1 = Point(sg.coords[0]), Point(sg.coords[-1])
                        if v.distance(e0) < MERGE_M or v.distance(e1) < MERGE_M:
                            continue                  # near its end vertex: the merge already decided
                        if best is None or d < best[0]:
                            best = (d, sg)
                    if best is not None:
                        q = best[1].interpolate(best[1].project(v))
                        table[c] = (q.x, q.y)
    if not table:
        return done
    for s, i in items:
        lab, g = done[s][i]
        parts = []
        for p in _polys(g):
            ext = [table.get((float(c[0]), float(c[1])), (float(c[0]), float(c[1]))) for c in p.exterior.coords]
            holes = [[table.get((float(c[0]), float(c[1])), (float(c[0]), float(c[1]))) for c in r.coords]
                     for r in p.interiors]
            q = _clean_poly(Polygon(ext, holes))
            if q is not None:
                parts.append(q)
        ng = _clean_poly(unary_union(parts) if len(parts) > 1 else (parts[0] if parts else None))
        if ng is not None and abs(ng.area - g.area) < 0.05 * max(g.area, 1.0):
            done[s][i] = (lab, ng)
    return done


PIN_M2 = 1.0            # empty ground this small inside the fabric is a pin-hole, not land
PIN_WIDTH_M = 0.3       # ...and so is any empty strip thinner than this


def _faces(done):
    """Rebuild every plot from ONE noded set of lines, so neighbours share their lines exactly.

    All plot outlines are noded together and polygonised into faces. Each face belongs to the
    plot it lies in (the larger one where two overlap by a hairline); an empty face that is a
    pin-hole goes to the plot bordering it most; real empty ground stays empty. A plot is the
    union of its faces, so two plots meet along the very same coordinates.
    """
    import shapely
    from shapely.ops import polygonize
    items = [(s, i) for s in done for i in range(len(done[s]))]
    geoms = [done[s][i][1] for s, i in items]
    tree = shapely.STRtree(geoms)
    lines = shapely.union_all([g.boundary for g in geoms], grid_size=GRID)      # snap-rounded noding
    owned = {k: [] for k in range(len(items))}
    pins = []
    for f in polygonize(lines):
        if f.area <= 1e-6:
            continue
        rp = f.representative_point()
        inside = [k for k in tree.query(rp) if geoms[k].contains(rp)]
        if inside:
            owned[max(inside, key=lambda k: geoms[k].area)].append(f)
        elif f.area < PIN_M2 or (f.length and 2.0 * f.area / f.length < PIN_WIDTH_M):
            pins.append(f)
    for f in pins:
        fb = f.buffer(0.01)
        cands = list(tree.query(fb))
        if not cands:
            continue
        k = max(cands, key=lambda k: geoms[k].boundary.intersection(fb).length)
        owned[k].append(f)
    for k, (s, i) in enumerate(items):
        if owned[k]:
            ng = _clean_poly(shapely.union_all(owned[k], grid_size=GRID))
            if ng is not None:
                done[s][i] = (done[s][i][0], ng)
    return done


# --------------------------------------------------------------------------- village driver

def tracker_status(village_code):
    z = zipfile.ZipFile(paths.TRACKER)
    sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8")
    shared = z.read("xl/sharedStrings.xml").decode("utf-8")
    z.close()
    vals = [re.sub(r"<.*?>", "", v) for v in re.findall(r"<si>(?:<t[^>]*>)?(.*?)(?:</t>)?</si>", shared, re.S)]
    vno = village_code.split("_")[-1].lstrip("0")
    out = {}
    for m in re.finditer(r'<row r="(\d+)".*?</row>', sheet, re.S):
        rid, body = m.group(1), m.group(0)
        g_ = re.search(r'<c r="G%s"[^>]*><v>(\d+)</v></c>' % rid, body)
        h_ = re.search(r'<c r="H%s"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        j_ = re.search(r'<c r="J%s"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        if g_ and h_ and j_ and g_.group(1) == vno:
            out[vals[int(h_.group(1))]] = vals[int(j_.group(1))]
    return out


def processing_order(village, bodies):
    """Railway land first (longest first), then outward ring by ring, larger first per ring."""
    rail = [s for s in bodies if s in engine._rail_parcels(village)]
    adj = {s: set() for s in bodies}
    for s in bodies:
        for n in neighbours.printed(village, s):
            if n in bodies and n != s:
                adj[s].add(n); adj[n].add(s)
    keys = list(bodies)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if bodies[keys[i]].distance(bodies[keys[j]]) <= 0.5:
                adj[keys[i]].add(keys[j]); adj[keys[j]].add(keys[i])
    order = sorted(rail, key=lambda s: -bodies[s].length)
    seen = set(order)
    frontier = list(order)
    while frontier:
        nxt = sorted({n for s in frontier for n in adj[s] if n not in seen}, key=lambda s: -bodies[s].area)
        order += nxt
        seen |= set(nxt)
        frontier = nxt
    order += sorted([s for s in bodies if s not in seen], key=lambda s: -bodies[s].area)
    return order, adj


def tidy_village(village, name, src, out):
    """Clean one village in place. Returns (gdf, report, order); writes `out` in EPSG:4326."""
    cur = gpd.read_file(src).to_crs(UTM)
    status = tracker_status(village)
    current = {}
    for s, sub in cur.groupby("survey_no"):
        cp = [(str(l or ""), _clean_poly(fmb_on_base._valid(g))) for l, g in zip(sub["subdiv_no"], sub.geometry)]
        current[str(s)] = [(l, g) for l, g in cp if g is not None]
    plots, report = {}, {}
    for s, cp in current.items():
        body = unary_union([g for _l, g in cp])
        outline = tidy_outline(body)
        interior, source = None, "your plots, tidied"
        pdf, drawn = regenerate(village, s, cp) if len(cp) > 1 else (None, None)
        if pdf:
            cand = partition(pdf, outline)
            if interior_gate(cand, cp):
                interior, source = cand, "your outline, PDF internal lines"
        if interior is None:
            interior = partition(cp, outline)
        plots[s] = absorb_slivers(interior)
        report[s] = {"source": source}
    bodies = {s: unary_union([g for _l, g in p]) for s, p in plots.items()}
    order, _adj = processing_order(village, bodies)
    done, finished = {}, []
    for s in order:
        near = [n for n in finished if bodies[n].distance(bodies[s]) <= 2.0 * SNAP_M]
        seniors = unary_union([unary_union([g for _l, g in done[n]]) for n in near]) if near else None
        p2, _moves = conflate(plots[s], seniors, tol=SNAP_M)
        done[s] = p2
        finished.append(s)
    report["_seams"] = close_seams(village, done)
    done = finish_coverage(done)
    rows = []
    for s in sorted(done, key=paths.survey_sort_key):
        mine = unary_union([g for _l, g in current[s]])
        new = unary_union([g for _l, g in done[s]])
        r = report[s]
        r["iou_yours"] = round(_iou(new, mine), 4)
        r["area_ratio"] = round(new.area / mine.area, 4) if mine.area else None
        r["max_move_m"] = round(_outer(new).hausdorff_distance(_outer(mine)), 2)
        m = _match(done[s], current[s])
        r["worst_plot"] = min(m, key=lambda x: x[1]) if m else None
        for lab, g in done[s]:
            g = set_precision(g, GRID)
            if g.is_empty:
                continue
            sub, kide = fmb_finish.survey_keys(s, lab)
            rows.append({"village_code": village, "village": name, "kide": kide, "survey_no": s,
                         "subdiv_no": sub, "status": status.get(s, "not in tracker"),
                         "shape_vs_yours": r["iou_yours"], "max_edge_move_m": r["max_move_m"],
                         "source": r["source"], "geometry": g})
    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs=UTM)
    gdf["area_sqm"] = gdf.geometry.area.round(1)
    gdf["area_acre"] = (gdf.geometry.area / ACRE).round(4)
    gdf = gdf[["village_code", "village", "kide", "survey_no", "subdiv_no", "area_acre", "area_sqm",
               "status", "shape_vs_yours", "max_edge_move_m", "source", "geometry"]]
    g4 = gdf.to_crs(4326)
    g4["geometry"] = [fmb_finish._write_safe(x) for x in g4.geometry]
    g4.to_file(out, driver="GeoJSON", COORDINATE_PRECISION=8)
    return gdf, report, order


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="village_clean", description=__doc__.splitlines()[0])
    ap.add_argument("village")
    ap.add_argument("name")
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    gdf, report, _order = tidy_village(a.village, a.name, a.src, a.out)
    print("%d plots, %d surveys -> %s" % (len(gdf), gdf.survey_no.nunique(), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
