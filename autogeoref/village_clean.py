"""Rebuild a hand-georeferenced village neatly: every survey from its PDF, one line per seam.

    python -m autogeoref.village_clean 35_04_074 Thirukatchur --src <geojson> --out <geojson>

Akash, 2026-09-24: internal lines were smudged and shapes had drifted after rounds of clip and
fill patching (vertex bloat, fragments, slivers, drawings off by up to 31 %). Patching again
cannot make it neat, so each survey is rebuilt:

  1. REGENERATE: every plot comes from the survey's PDF sheet (straight internal lines, exact
     legal shape, the sheet's plot count), placed rigidly where the survey currently sits. The
     pose maximises agreement with the current outline AND the current plots by number, so a
     rectangle cannot come back turned 180 degrees with its plots mirrored.
  2. CONFLATE, railway land first then outward: each survey is fitted to the already-finished
     surveys it borders (printed PDF neighbours, or ones it touches or overlaps). Only OUTER
     boundary points move: a point within SNAP_TOL_M of a finished neighbour's boundary moves
     onto it (onto a vertex when one is close, so the two share it); then every edge whose two
     ends sit on that boundary is replaced by the neighbour's own boundary between them, so the
     seam is one line, not two nearly-parallel ones. Internal points never move, so internal
     lines stay the PDF's straight lines. Anything still overlapping is cut along the neighbour's
     edge; fragments under FRAGMENT_M2 are dropped.
  3. The report says, per survey, how far any boundary point moved and how close the result
     stays to the drawing.
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


def rebuild_village(village, name, src, out):
    cur = gpd.read_file(src).to_crs(UTM)
    status = tracker_status(village)
    current = {s: [(str(l or ""), fmb_on_base._valid(g)) for l, g in zip(sub["subdiv_no"], sub.geometry)]
               for s, sub in cur.groupby("survey_no")}
    regen, drawn, report = {}, {}, {}
    for s, cp in current.items():
        plots, outline = regenerate(village, s, cp)
        if plots is None:
            regen[s] = cp                      # no sheet: keep as placed
            report[s] = {"note": "no PDF sheet - kept as placed"}
        else:
            regen[s], drawn[s] = plots, outline
    bodies = {s: unary_union([g for _l, g in p]) for s, p in regen.items()}
    order, adj = processing_order(village, bodies)
    done, finished = {}, []
    for s in order:
        senior_names = [n for n in finished if n in adj[s] or bodies[n].intersects(bodies[s])]
        seniors = unary_union([unary_union([g for _l, g in done[n]]) for n in senior_names]) if senior_names else None
        plots, moves = conflate(regen[s], seniors)
        # a seam may only close over EMPTY ground: whatever the survey gained that another survey's
        # drawing occupies is handed back, even if that survey is not processed yet (53 closed its
        # rail seam straight across 47B's 2.6 m strip on 2026-09-24)
        pending = [n for n in bodies if n != s and n not in done and bodies[n].distance(bodies[s]) <= SNAP_TOL_M]
        if pending:
            own = unary_union([g for _l, g in regen[s]])
            reserved = unary_union([bodies[n] for n in pending]).difference(own)
            if not reserved.is_empty:
                plots = [(l, _clean_poly(g.difference(reserved)) if g.intersects(reserved) else g) for l, g in plots]
                plots = [(l, g) for l, g in plots if g is not None]
        # anything overlapping a finished survey that is NOT a neighbour is still cut
        others = [n for n in finished if n not in senior_names]
        if others:
            ob = unary_union([unary_union([g for _l, g in done[n]]) for n in others])
            plots = [(l, _clean_poly(g.difference(ob)) if g.intersects(ob) else g) for l, g in plots]
            plots = [(l, g) for l, g in plots if g is not None]
        # when the cut costs this survey too much, settle each overlap by proportion instead: the
        # side that loses the SMALLER share of itself gives the ground up (53 gives 8 % so that
        # 47B does not lose 92 % of itself; 2026-09-24)
        drawn_area = sum(g.area for _l, g in regen[s])
        kept_area = sum(g.area for _l, g in plots)
        if drawn_area and (1.0 - kept_area / drawn_area > MAX_LOSS or len(plots) < len(regen[s])):
            mine = unary_union([g for _l, g in regen[s]])
            keep_out = []
            yielded = []
            for n in finished:
                nb = unary_union([g for _l, g in done[n]])
                ov = nb.intersection(mine).area
                if ov <= 0.05:
                    continue
                if ov / nb.area < ov / mine.area and ov / nb.area <= MAX_LOSS:
                    done[n] = [(l, _clean_poly(g.difference(mine)) if g.intersects(mine) else g)
                               for l, g in done[n]]
                    done[n] = [(l, g) for l, g in done[n] if g is not None]
                    yielded.append(n)
                else:
                    keep_out.append(nb)
            plots = regen[s]
            if keep_out:
                ko = unary_union(keep_out)
                plots = [(l, _clean_poly(g.difference(ko)) if g.intersects(ko) else g) for l, g in plots]
                plots = [(l, g) for l, g in plots if g is not None]
            if yielded:
                report.setdefault(s, {})["note"] = "%s yielded the overlap" % ",".join(yielded)
        done[s] = plots
        finished.append(s)
        body = unary_union([g for _l, g in plots])
        r = report.setdefault(s, {})
        r["max_move_m"] = round(max(moves), 2) if moves else 0.0
        r["iou_pdf"] = round(_iou(body, drawn[s]), 3) if s in drawn else None
    rows = []
    for s in sorted(done, key=paths.survey_sort_key):
        for lab, g in done[s]:
            g = set_precision(g, GRID)
            if g.is_empty:
                continue
            sub, kide = fmb_finish.survey_keys(s, lab)
            rows.append({"village_code": village, "village": name, "kide": kide, "survey_no": s,
                         "subdiv_no": sub, "status": status.get(s, "not in tracker"),
                         "shape_iou_pdf": report[s].get("iou_pdf"),
                         "max_edge_move_m": report[s].get("max_move_m"),
                         "source": "PDF sheet, rigid, seams conflated" if s in drawn else "as placed",
                         "geometry": g})
    gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs=UTM)
    gdf["area_sqm"] = gdf.geometry.area.round(1)
    gdf["area_acre"] = (gdf.geometry.area / ACRE).round(4)
    gdf = gdf[["village_code", "village", "kide", "survey_no", "subdiv_no", "area_acre", "area_sqm",
               "status", "shape_iou_pdf", "max_edge_move_m", "source", "geometry"]]
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
    gdf, report, _order = rebuild_village(a.village, a.name, a.src, a.out)
    print("%d plots, %d surveys -> %s" % (len(gdf), gdf.survey_no.nunique(), a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
