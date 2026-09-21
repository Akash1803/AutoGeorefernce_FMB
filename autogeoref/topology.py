"""Topology between parcels the tool may move. Anchors are never touched.

No snapping. The 2026-09-18 16:20 run used a 0.30 m snap across every polygon of the village; on
526A, whose 177 subdivision plots include slivers under a metre across, that collapsed two plots to
nothing and cut 153 m2 off a third. Errors are fixed only by clipping an overlap out of one side and
by giving a gap to the polygon that borders it most, and every edit is refused if it would empty a
polygon or take more than GUARD_FRACTION of its area.
"""
import csv
import itertools
import math

from shapely.geometry import LineString, Point, Polygon
from shapely.ops import nearest_points, unary_union

from . import paths

MIN_AREA = 0.02          # square metres: smaller than this is numerical noise
GUARD_FRACTION = 0.25    # refuse any edit that takes more than this share of a polygon


MITRE = dict(join_style=2, mitre_limit=5.0)
CLEAN_TOL = 0.02         # metres: hairline gaps closed, near-collinear vertices dropped
MIN_WIDTH = 0.20         # metres: anything thinner than this is a sliver, not a plot

CONFORM_TOL = 0.60       # metres: two automated parcels this close share a boundary and must draw it
                         # with identical coordinates (adjustment residuals are 0.2 m rms, 0.4 m max)
ANCHOR_CONFORM_TOL = 1.50  # against the team's own parcel, which is affine: 171 sat 1.25 m off 48A
ON_LINE = 0.02           # metres: a vertex this close to a boundary is on it
INSERT_DEG = 20.0        # a plot edge within this of parallel to a settled boundary runs along it
MERGE_TOL = 0.05         # metres: two plot vertices this close are the same drawn point
_KEY = 6                 # decimals for coordinate keys


def _mitre(geom, d):
    return geom.buffer(d, **MITRE)


def clean(geom, tol=CLEAN_TOL, min_width=MIN_WIDTH):
    """Tidy a topology-edited polygon so it looks like a drawn parcel, not a computed one.

    Closes hairline gaps between a filled strip and its plot (they otherwise draw as a double
    line), removes needles thinner than `min_width`, keeps the largest piece, and drops the
    near-collinear vertices that clipping leaves behind. Mitred buffers keep every real corner
    sharp; opening a polygon by `min_width/2` and closing it again returns the original wherever
    it is wider than `min_width`.
    """
    if geom is None or geom.is_empty:
        return geom
    g = geom.buffer(0)
    g = _mitre(_mitre(g, tol), -tol)                          # close hairlines
    opened = _mitre(_mitre(g, -min_width / 2.0), min_width / 2.0)   # remove needles
    if not opened.is_empty and opened.area >= 0.98 * g.area:
        g = opened
    if g.geom_type == "MultiPolygon":
        g = max(g.geoms, key=lambda p: p.area)
    elif g.geom_type == "GeometryCollection":
        polys = [p for p in g.geoms if p.geom_type in ("Polygon", "MultiPolygon")]
        g = max(unary_union(polys).geoms, key=lambda p: p.area) if polys and unary_union(polys).geom_type == "MultiPolygon" else (unary_union(polys) if polys else Polygon())
    g = g.simplify(tol, preserve_topology=True)
    return g.buffer(0)


def _largest(geom):
    """Repair a geometry without silently discarding parts of a multipolygon."""
    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == "GeometryCollection":
        polys = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
        return unary_union(polys) if polys else Polygon()
    return geom


def check(geoms):
    keys = sorted(geoms, key=paths.survey_sort_key)
    overlaps = []
    for a, b in itertools.combinations(keys, 2):
        area = geoms[a].intersection(geoms[b]).area
        if area > MIN_AREA:
            overlaps.append([a, b, round(area, 2)])
    union = unary_union([g.buffer(0) for g in geoms.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    gaps = [round(Polygon(r).area, 2) for p in parts for r in p.interiors if Polygon(r).area > 0.01]
    return {"overlap_pairs": len(overlaps), "overlap_sqm": round(sum(o[2] for o in overlaps), 2),
            "overlaps": overlaps, "gaps": len(gaps), "gap_sqm": round(sum(gaps), 2),
            "gap_list": sorted(gaps, reverse=True)[:15],
            "invalid": [k for k, g in geoms.items() if not g.is_valid],
            "union_parts": len(parts)}


def fix(geoms, movable, rail, gap_tol=1.20):
    """Clip overlaps and fill gaps between `movable` parcels. Returns (geoms, report)."""
    out = {k: g.buffer(0) for k, g in geoms.items()}
    report = {"clips": [], "fills": [], "refused": [], "rail_conflicts": []}
    keys = sorted(out, key=paths.survey_sort_key)

    def edit(key, new_geom, why, taken):
        old = out[key]
        if new_geom is None or new_geom.is_empty:
            report["refused"].append([key, why, round(taken, 2), "would empty the polygon"])
            return False
        if old.area > 0 and taken / old.area > GUARD_FRACTION:
            report["refused"].append([key, why, round(taken, 2),
                                      "would take %.0f %% of a %.2f m2 polygon"
                                      % (100 * taken / old.area, old.area)])
            return False
        out[key] = _largest(new_geom)
        return True

    # 1. overlaps: one side yields
    for _round in range(4):
        changed = 0
        for a, b in itertools.combinations(keys, 2):
            inter = out[a].intersection(out[b])
            if inter.is_empty or inter.area <= MIN_AREA:
                continue
            ra, rb = a in rail, b in rail
            if ra != rb:
                # railway land was carved out of the older parcel, so the strip keeps its shape.
                # A parcel the tool placed is clipped to it (guarded); a team parcel is only reported.
                pair = (a, b, round(inter.area, 2))
                if pair not in report["rail_conflicts"]:
                    report["rail_conflicts"].append(pair)
                other = b if ra else a
                strip = a if ra else b
                if other in movable and edit(other, out[other].difference(out[strip]),
                                             "clip against railway land %s" % strip, inter.area):
                    report["clips"].append([other, "clipped against railway land", strip, round(inter.area, 2)])
                    changed += 1
                continue
            big, small = (a, b) if out[a].area >= out[b].area else (b, a)
            if big not in movable:
                big, small = small, big
            if big not in movable:
                continue
            if edit(big, out[big].difference(out[small]), "clip against %s" % small, inter.area):
                report["clips"].append([big, "clipped against", small, round(inter.area, 2)])
                changed += 1
        if not changed:
            break

    # 2. enclosed gaps, then thin open gaps between two surveys
    union = unary_union([g.buffer(0) for g in out.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    candidates = [(Polygon(r), "enclosed gap") for p in parts for r in p.interiors
                  if Polygon(r).area > MIN_AREA]
    for a, b in itertools.combinations(keys, 2):
        d = out[a].distance(out[b])
        if d <= 0 or d > gap_tol:
            continue
        # mitred buffers: round ones left arcs on the filled pieces (the curve Akash saw on 2026-09-20)
        strip = (_mitre(out[a], gap_tol).intersection(_mitre(out[b], gap_tol))
                 .difference(out[a]).difference(out[b]))
        # only the part between the two parcels: the buffers also overlap in square "ears" past
        # the ends of the gap, and those ears are not a gap
        strip = strip.intersection(unary_union([out[a], out[b]]).convex_hull)
        for piece in (strip.geoms if strip.geom_type.startswith("Multi") else [strip]):
            if piece.is_empty or piece.area <= MIN_AREA or piece.length == 0:
                continue
            if 2 * piece.area / piece.length > gap_tol:       # mean width: must be a sliver
                continue
            candidates.append((piece, "gap between %s and %s" % (a, b)))
    for gap, kind in candidates:
        best, best_share = None, 0.0
        for key in movable:
            share = out[key].buffer(0.05).intersection(gap).area
            if share > best_share:
                best, best_share = key, share
        if best is None:
            continue
        out[best] = _largest(unary_union([out[best], gap]))
        report["fills"].append([round(gap.area, 2), "into", best, kind])
    edited = {r[0] for r in report["clips"]} | {r[2] for r in report["fills"]}
    for key in edited:
        if key in movable:
            out[key] = clean(out[key])        # tidy only what was edited; anchors are never here
    return out, report


def apply_to_parts(parts, fixed_body, original_body):
    """Carry a survey's topology-fixed body down to its sub-plots.

    `parts` is [(props, geom)] as placed. Every part is clipped to the fixed body (so a clipped
    overlap disappears from the plot that held it) and each filled gap piece goes to the part
    that borders it most. Returns [(props, geom, note)] with note "" / "clipped" / "filled".
    """
    out = []
    added = fixed_body.difference(original_body)
    pieces = [p for p in (added.geoms if added.geom_type.startswith("Multi") else [added])
              if not p.is_empty and p.area > MIN_AREA]
    owner = {}
    for i, piece in enumerate(pieces):
        best, best_len = None, 0.0
        for k, (_props, g) in enumerate(parts):
            shared = g.buffer(0.05).intersection(piece).area
            if shared > best_len:
                best, best_len = k, shared
        if best is not None:
            owner.setdefault(best, []).append(piece)
    for k, (props, g) in enumerate(parts):
        ng = g.intersection(fixed_body)
        note = ""
        if abs(ng.area - g.area) > MIN_AREA:
            note = "clipped"
        if k in owner:
            ng = unary_union([ng] + owner[k])
            note = (note + ",filled").strip(",")
        ng = _largest(ng)
        if note:
            ng = clean(ng)                 # untouched plots stay exactly as the sheet drew them
        if ng is None or ng.is_empty or ng.geom_type != "Polygon":
            ng, note = g, "kept (edit would empty the plot)"
        out.append((props, ng, note))
    return out


def _key(c):
    return (round(c[0], _KEY), round(c[1], _KEY))


def _boundary_lines(body):
    polys = body.geoms if body.geom_type == "MultiPolygon" else [body]
    lines = []
    for p in polys:
        lines.append(LineString(p.exterior.coords))
        for r in p.interiors:
            lines.append(LineString(r.coords))
    return lines


def _outer(geoms):
    body = unary_union([g.buffer(0) for g in geoms])
    if body.geom_type == "MultiPolygon":
        body = max(body.geoms, key=lambda p: p.area)
    elif body.geom_type == "GeometryCollection":
        body = _largest(body)
    return body, list(body.exterior.coords)[:-1]


def _dedupe(coords):
    out = []
    for c in coords:
        if not out or _key(c) != _key(out[-1]):
            out.append(c)
    if len(out) > 1 and _key(out[0]) == _key(out[-1]):
        out.pop()
    return out


def _with_ring(geom, exterior, interiors):
    try:
        return Polygon(exterior, interiors)
    except (ValueError, TypeError):
        return None


def node_parts(geoms, merge_tol=MERGE_TOL, on_line=ON_LINE):
    """Make a parcel's plots share every vertex they touch.

    Converted sheets draw T-junctions (a plot corner on a neighbour plot's edge with no vertex
    there) and near-duplicate corners a few centimetres apart. Both are harmless until a vertex is
    moved: on 2026-09-21 conform moved a 46A corner 0.15 m and the neighbouring plot's straight edge
    stayed, opening a 100 m hairline between two plots of the same parcel. Here vertices within
    `merge_tol` become one point and a vertex lying on another plot's edge is inserted into it.
    Returns (new_geoms, changed).
    """
    canon = {}
    points = []
    for g in geoms:
        for ring in [g.exterior, *g.interiors]:
            for c in list(ring.coords)[:-1]:
                k = _key(c)
                if k in canon:
                    continue
                hit = next((pc for pc in points if math.hypot(pc[0] - c[0], pc[1] - c[1]) <= merge_tol), None)
                canon[k] = hit if hit is not None else c
                if hit is None:
                    points.append(c)

    def merged(coords):
        return _dedupe([canon.get(_key(c), c) for c in coords])

    rings = [[merged(list(r.coords)[:-1]) for r in [g.exterior, *g.interiors]] for g in geoms]
    all_pts = [Point(c) for c in points]
    out, changed = [], False
    for gi, g in enumerate(geoms):
        new_rings = []
        for cs in rings[gi]:
            res = []
            m = len(cs)
            for i in range(m):
                a, b = cs[i], cs[(i + 1) % m]
                res.append(a)
                seg = LineString([a, b])
                L = seg.length
                if L <= on_line:
                    continue
                ins = []
                for q in all_pts:
                    if seg.distance(q) > on_line:
                        continue
                    prm = seg.project(q) / L
                    if 1e-6 < prm < 1 - 1e-6 and _key((q.x, q.y)) not in (_key(a), _key(b)):
                        ins.append((prm, (q.x, q.y)))
                ins.sort()
                res.extend(c for _pm, c in ins)
            new_rings.append(_dedupe(res))
        ext, ints = new_rings[0], [r for r in new_rings[1:] if len(r) >= 3]
        ng = Polygon(ext, ints) if len(ext) >= 3 else g
        if not ng.is_valid:
            ng = g
        if not ng.equals_exact(g, 1e-9):
            changed = True
        out.append(ng)
    return out, changed


def conform(parts_by_survey, settled, order, anchors=(), tol=CONFORM_TOL, anchor_tol=ANCHOR_CONFORM_TOL):
    """Make every shared boundary identical, vertex for vertex.

    Clipping and filling with buffers leaves two boundaries a centimetre apart: on 2026-09-21 every
    pair of neighbouring parcels of the ring-2 result crossed each other one to seven times (the
    X-shaped double lines Akash pointed at). Here each movable parcel, most certain first, is
    conformed to the parcels already settled around it, plot by plot:

    - a plot vertex within tolerance of a settled boundary moves onto it; a settled corner within
      reach wins over the nearer edge (at 48A's corner 47B was pulled onto 171's line, 0.3 m away,
      instead of onto 48A's corner 1.2 m away, and a 0.4 m strip stayed open along the anchor);
    - a settled vertex lying beside a plot edge that runs along that boundary (within INSERT_DEG of
      parallel) is inserted into the edge, so both rings carry the same points; the edge's own
      endpoints need not be on any line (171 had no vertex near 48A's corner at all).

    A vertex shared by two plots moves identically in both, and an inserted point goes into every
    ring holding that edge, so a parcel's plots stay stitched. Nothing settled is ever changed; a
    plot whose edit would be invalid or would change its area by more than GUARD_FRACTION keeps its
    geometry and is reported.

    parts_by_survey: {survey: [(props, geom)]} for the movable parcels.
    settled: {survey: body} for the fixed parcels (the team's own geometry).
    order: movable surveys, most certain first.
    Returns ({survey: [(props, geom, note)]}, report).
    """
    settled = dict(settled)
    out, report = {}, {"conformed": {}, "refused": []}
    cos_lim = math.cos(math.radians(INSERT_DEG))
    for s in order:
        originals = parts_by_survey[s]
        noded, _was_noded = node_parts([g for _p, g in originals])
        parts = [(p, g) for (p, _g0), g in zip(originals, noded)]
        geoms = [g for _p, g in parts]
        body = _body_of(geoms)
        refs = []
        for n, nb in settled.items():
            if n == s or nb is None or nb.is_empty:
                continue
            reach = anchor_tol if n in anchors else tol
            if nb.distance(body) > reach:
                continue
            lines = _boundary_lines(nb)
            vdirs = {}
            for l in lines:
                cs = list(l.coords)[:-1]
                m = len(cs)
                for i, c in enumerate(cs):
                    a, b = cs[i - 1], cs[(i + 1) % m]
                    dirs = []
                    for o in (a, b):
                        dx, dy = o[0] - c[0], o[1] - c[1]
                        L = math.hypot(dx, dy)
                        if L > 1e-9:
                            dirs.append((dx / L, dy / L))
                    vdirs.setdefault(_key(c), (Point(c), []))[1].extend(dirs)
            refs.append((unary_union(lines), vdirs, reach, n, n in anchors))
        if not refs:
            out[s] = [(p, g, "" if g.equals_exact(g0, 1e-9) else "noded")
                      for (p, g), (_p0, g0) in zip(parts, originals)]
            settled[s] = body
            continue

        # every vertex and edge of every plot
        keys, edges = {}, {}
        for g in geoms:
            for ring in [g.exterior, *g.interiors]:
                cs = list(ring.coords)[:-1]
                for i, c in enumerate(cs):
                    keys.setdefault(_key(c), c)
                    ka, kb = _key(c), _key(cs[(i + 1) % len(cs)])
                    edges[(ka, kb) if ka <= kb else (kb, ka)] = (c, cs[(i + 1) % len(cs)])

        # 1. where each vertex goes: a settled corner beats a settled edge, an anchor beats a peer
        vmap, how = {}, {}
        for k, c in keys.items():
            pt = Point(c)
            cands = []
            for lines, vdirs, reach, n, is_anchor in refs:
                d = lines.distance(pt)
                if d > reach:
                    continue
                near = [(pt.distance(q), q) for q, _d in vdirs.values() if pt.distance(q) <= reach]
                if near:
                    dv, q = min(near, key=lambda x: x[0])
                    cands.append((0, 0 if is_anchor else 1, dv, (q.x, q.y), n))
                else:
                    q = nearest_points(lines, pt)[0]
                    cands.append((1, 0 if is_anchor else 1, d, (q.x, q.y), n))
            if not cands:
                continue
            best = min(cands)
            if best[2] <= 1e-9:
                continue
            vmap[k], how[k] = best[3], best
        # a settled corner takes one plot vertex only, the nearest; the rest go onto the edge
        # (46A's plot corners along 47B all went to one 47B vertex and were then dropped, 2026-09-21)
        by_target = {}
        for k, h in how.items():
            if h[0] == 0:
                by_target.setdefault(_key(h[3]), []).append((h[2], k))
        for tk, lst in by_target.items():
            lst.sort()
            for _d, k in lst[1:]:
                lines = next(r[0] for r in refs if r[3] == how[k][4])
                q = nearest_points(lines, Point(keys[k]))[0]
                vmap[k] = (q.x, q.y)
                how[k] = (1, how[k][1], lines.distance(Point(keys[k])), (q.x, q.y), how[k][4])
        # two ends of a real edge must not land on the same point
        for (ka, kb), (a, b) in edges.items():
            if ka in vmap and kb in vmap and _key(vmap[ka]) == _key(vmap[kb]) \
                    and Point(a).distance(Point(b)) > tol:
                drop = ka if how[ka][2] >= how[kb][2] else kb
                del vmap[drop], how[drop]

        def moved(coords):
            return _dedupe([vmap.get(_key(c), c) for c in coords])

        # 2. settled vertices beside an edge that runs along their boundary are inserted into it
        new_keys = {_key(vmap.get(k, c)) for k, c in keys.items()}
        inserts = {}                                       # (key_u, key_w) after moving -> [(param, coord)]
        seg_of = {}
        for _ek, (a, b) in edges.items():
            # (a, b) is the ring's own order; the dict key is sorted and must not be paired with it
            u, w = vmap.get(_key(a), a), vmap.get(_key(b), b)
            if _key(u) == _key(w):
                continue
            seg = LineString([u, w])
            ku, kw = _key(u), _key(w)
            # the segment keeps its own direction (u -> w); insertion parameters refer to it
            seg_of[(ku, kw) if ku <= kw else (kw, ku)] = (seg, ku, kw)
        for lines, vdirs, reach, n, _is_anchor in refs:
            for rk, (q, dirs) in vdirs.items():
                if rk in new_keys or not dirs:
                    continue
                best = None
                for ek, (seg, _ku, _kw) in seg_of.items():
                    d = seg.distance(q)
                    if d > reach:
                        continue
                    (x0, y0), (x1, y1) = seg.coords[0], seg.coords[-1]
                    L = seg.length
                    ex, ey = (x1 - x0) / L, (y1 - y0) / L
                    if max(abs(ex * dx + ey * dy) for dx, dy in dirs) < cos_lim:
                        continue                           # the edge does not run along this boundary
                    param = seg.project(q) / L
                    if param <= 1e-6 or param >= 1 - 1e-6:
                        continue
                    if best is None or d < best[0]:
                        best = (d, ek, param)
                if best is not None:
                    inserts.setdefault(best[1], []).append((best[2], (q.x, q.y)))
        for ek in inserts:
            inserts[ek].sort()

        def inserted(coords):
            res = []
            m = len(coords)
            for i in range(m):
                a, b = coords[i], coords[(i + 1) % m]
                res.append(a)
                ka, kb = _key(a), _key(b)
                ek = (ka, kb) if ka <= kb else (kb, ka)
                if ek in inserts:
                    _seg, ku, kw = seg_of[ek]
                    forward = (ka == ku and kb == kw)          # this ring walks the segment u -> w
                    pts = inserts[ek] if forward else list(reversed(inserts[ek]))
                    res.extend(c for _pm, c in pts)
            return _dedupe(res)

        # 3. apply to every plot; validate each
        new_parts, n_moved, max_move, n_ins = [], 0, 0.0, 0
        for (props, g), (_p0, g_orig) in zip(parts, originals):
            ext0 = list(g.exterior.coords)[:-1]
            ext = inserted(moved(ext0))
            ints = [inserted(moved(list(r.coords)[:-1])) for r in g.interiors]
            ng = _with_ring(g, ext, [r for r in ints if len(r) >= 3]) if len(ext) >= 3 else None
            if ng is not None and not ng.is_valid:
                ng = ng.buffer(0)
                if ng.geom_type == "MultiPolygon":
                    ng = max(ng.geoms, key=lambda p: p.area)
                elif ng.geom_type == "GeometryCollection":
                    ng = _largest(ng)
            changed = ng is not None and not ng.equals_exact(g, 1e-9)
            if ng is None or ng.is_empty or ng.geom_type != "Polygon" or \
                    (g.area > 0 and abs(ng.area - g.area) / g.area > GUARD_FRACTION):
                report["refused"].append([s, props.get("poly_id"), "plot edit invalid or too large"])
                new_parts.append((props, g_orig, "kept (conform refused)"))
                continue
            if changed:
                mv = [Point(c).distance(Point(vmap[_key(c)])) for c in ext0 if _key(c) in vmap]
                n_moved += len(mv)
                max_move = max([max_move] + mv)
                n_ins += max(0, len(ext) - len(_dedupe(moved(ext0))))
            note = "conformed" if changed else ("" if ng.equals_exact(g_orig, 1e-9) else "noded")
            new_parts.append((props, ng, note))
        out[s] = new_parts
        settled[s] = _body_of([g for _p, g, _n in new_parts])
        if n_moved or n_ins:
            report["conformed"][s] = {"vertices_moved": n_moved, "max_move_m": round(max_move, 3),
                                      "vertices_inserted": n_ins}
    return out, report

MAX_FILL_M2 = 25.0       # an enclosed gap larger than this is a missing parcel, not a sliver


def _body_of(parts):
    return unary_union([g.buffer(0) for g in parts if g is not None and not g.is_empty])


def resolve(parts_by_survey, settled, order, rail=(), anchors=None, tol=CONFORM_TOL,
            anchor_tol=ANCHOR_CONFORM_TOL):
    """The whole topology stage at plot level, with exact geometry only.

    1. conform: shared boundaries become identical, vertex for vertex (see `conform`).
    2. clip: an overlap that is still there is a real crossing beyond tolerance; the less certain
       parcel yields, plot by plot, by exact difference. Railway land keeps its shape and the
       ordinary parcel yields to it; an overlap with the team's parcel is taken out of ours.
    3. fill: a hole enclosed by the parcels that is small enough to be a sliver goes, by exact
       union, to the plot that borders it most.
    No buffers anywhere: the 2026-09-21 ring-2 files showed that buffered clips, fills and
    cleaning leave every pair of neighbours crossing each other, and that per-plot cleaning
    breaks the vertices adjacent plots share.

    Returns ({survey: [(props, geom, note)]}, report).
    """
    anchors = set(anchors if anchors is not None else settled)
    out, report = conform(parts_by_survey, settled, order, anchors=anchors, tol=tol, anchor_tol=anchor_tol)
    report.update({"clips": [], "fills": [], "rail_conflicts": []})
    rank = {s: i for i, s in enumerate(order)}
    bodies = dict(settled)
    bodies.update({s: _body_of([g for _p, g, _n in out[s]]) for s in out})

    def clip_survey(s, against, why):
        taken_total = 0.0
        new = []
        for props, g, note in out[s]:
            inter = g.intersection(bodies[against])
            if inter.is_empty or inter.area <= MIN_AREA:
                new.append((props, g, note))
                continue
            ng = _largest(g.difference(bodies[against]))
            if ng is None or ng.is_empty or ng.geom_type != "Polygon" or inter.area / g.area > GUARD_FRACTION:
                report["refused"].append([s, props.get("poly_id"), why, round(inter.area, 2)])
                new.append((props, g, note))
                continue
            taken_total += inter.area
            new.append((props, ng, ",".join(x for x in (note, "clipped") if x)))
        if taken_total > 0:
            out[s] = new
            bodies[s] = _body_of([g for _p, g, _n in out[s]])
            report["clips"].append([s, why, against, round(taken_total, 2)])
        return taken_total

    # 2. residual overlaps
    for _round in range(3):
        changed = 0.0
        keys = sorted(bodies, key=lambda k: (rank.get(k, -1), paths.survey_sort_key(k)))
        for a, b in itertools.combinations(keys, 2):
            inter = bodies[a].intersection(bodies[b]).area
            if inter <= MIN_AREA:
                continue
            ma, mb = a in out, b in out
            if not ma and not mb:
                continue                                  # two team parcels: theirs to settle
            ra, rb = a in rail, b in rail
            if ra != rb:
                pair = (a, b, round(inter, 2))
                if pair not in report["rail_conflicts"]:
                    report["rail_conflicts"].append(pair)
                other, strip = (b, a) if ra else (a, b)
                if other in out:
                    changed += clip_survey(other, strip, "clip against railway land %s" % strip)
                continue
            if ma and mb:
                yielder, keeper = (b, a) if rank[a] <= rank[b] else (a, b)   # the less certain yields
            else:
                yielder, keeper = (a, b) if ma else (b, a)                  # ours yields to the team's
            changed += clip_survey(yielder, keeper, "clip against %s" % keeper)
        if changed <= MIN_AREA:
            break

    # 3. enclosed slivers
    union = unary_union(list(bodies.values()))
    polys = union.geoms if union.geom_type == "MultiPolygon" else [union]
    for p in polys:
        for r in p.interiors:
            hole = Polygon(r)
            if hole.area <= MIN_AREA or hole.area > MAX_FILL_M2:
                continue
            best, best_len = None, 0.0
            hb = hole.boundary
            for s in out:
                for k, (props, g, note) in enumerate(out[s]):
                    # the union nodes coordinates at 1e-10, so an exact line intersection can come
                    # back as points; measure the shared length within a hair instead
                    shared = g.boundary.buffer(0.01, cap_style=2).intersection(hb).length
                    if shared > best_len:
                        best, best_len = (s, k), shared
            if best is None or best_len <= 0:
                continue
            s, k = best
            props, g, note = out[s][k]
            ng = unary_union([g, hole])
            if ng.geom_type != "Polygon":
                ng = _largest(ng)
            if ng is None or ng.geom_type != "Polygon":
                continue
            out[s][k] = (props, ng, ",".join(x for x in (note, "filled") if x))
            bodies[s] = _body_of([g2 for _p, g2, _n in out[s]])
            report["fills"].append([round(hole.area, 2), "into", s, props.get("poly_id")])
    return out, report


MAX_LOSS = 0.05                 # a parcel that lost more than this of its rigid area was mis-placed
SPIKE_DEG = 10.0
SPIKE_ARM_M = 2.0
SIBLING_OVERLAP_M2 = 0.5


def _spiky(geom):
    ring = list(geom.exterior.coords)
    n = len(ring) - 1
    for i in range(n):
        a, b, c = ring[i - 1], ring[i], ring[(i + 1) % n]
        v1 = (a[0] - b[0], a[1] - b[1])
        v2 = (c[0] - b[0], c[1] - b[1])
        l1, l2 = math.hypot(*v1), math.hypot(*v2)
        if min(l1, l2) < SPIKE_ARM_M:
            continue
        ang = math.degrees(math.acos(max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1]) / (l1 * l2)))))
        if ang < SPIKE_DEG:
            return True
    return False


def conform_as_body(survey, rigid_parts, settled, anchors=(), tol=CONFORM_TOL, anchor_tol=ANCHOR_CONFORM_TOL):
    """Conform one survey as a single body, then carry the new outline down to its plots.

    Used when plot-by-plot conform pulled a survey's plots apart. The dissolved outline is
    conformed against `settled` (anchors and earlier outputs), each plot is clipped to it and
    every filled piece goes to the plot bordering it most. Returns (parts, note) or (None, why).
    """
    body = unary_union([g for _p, g in rigid_parts])
    if body.is_empty or body.geom_type != "Polygon":
        return None, "body not a single polygon"
    out, rep = conform({survey: [({"poly_id": 0, "plot_no": survey}, body)]}, settled, [survey],
                       anchors=set(anchors), tol=tol, anchor_tol=anchor_tol)
    new_body = out.get(survey, [(None, None, None)])[0][1]
    if new_body is None or new_body.is_empty or new_body.geom_type != "Polygon" or not new_body.is_valid:
        return None, "body conform failed"
    for s, g in settled.items():
        if new_body.intersection(g).area > SIBLING_OVERLAP_M2:
            return None, "body conform overlaps %s" % s
    parts = apply_to_parts(rigid_parts, new_body, body)
    moved = rep.get("conformed", {}).get(survey, {}).get("max_move_m", 0.0)
    return [(p, g, ("body-conformed," + n).strip(",")) for p, g, n in parts], "conformed as one body, max %.2f m" % moved


def sanity(rigid_parts, new_parts, max_loss=MAX_LOSS):
    """Guard the resolved plots of one survey before they are written.

    Returns (parts, notes). A plot that came out invalid, spiky (a corner under SPIKE_DEG with
    arms over SPIKE_ARM_M) or overlapping a sibling plot is replaced by its rigid geometry clipped
    to the resolved body ("rigid-clipped"). If the survey lost more than `max_loss` of its rigid
    area to clipping, the note "topology cut N %" is added: the placement, not the topology, is
    at fault, and the caller marks the parcel red. 51 lost 36 % and 50A 7 % on 2026-09-21 and
    were delivered with sliver plots.
    """
    rigid = {i: g for i, (_p, g) in enumerate(rigid_parts)}
    body_rigid = unary_union([g for _p, g in rigid_parts])
    body_new = unary_union([g for _p, g, _n in new_parts])
    notes = []
    loss = (body_rigid.area - body_new.area) / body_rigid.area if body_rigid.area else 0.0
    if loss > max_loss:
        notes.append("topology cut %.0f %%" % (100 * loss))
    if body_new.geom_type != "Polygon" and body_rigid.geom_type == "Polygon":
        # the plots came apart (one conformed, one refused): the sheet is one piece, keep it so
        return [(props, rigid[i], "rigid (plots came apart)") for i, (props, _g, _n) in enumerate(new_parts)], \
            notes + ["plots came apart, survey kept rigid"]
    parts = list(new_parts)
    bad = set()
    for i, (_p, g, _n) in enumerate(parts):
        if g is None or g.is_empty or g.geom_type != "Polygon" or not g.is_valid or _spiky(g):
            bad.add(i)
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            gi, gj = parts[i][1], parts[j][1]
            if gi is not None and gj is not None and gi.intersection(gj).area > SIBLING_OVERLAP_M2:
                bad.add(i if gi.area < gj.area else j)
    for i in sorted(bad):
        props, _g, note = parts[i]
        fallback = _largest(rigid[i].intersection(body_new)) if i in rigid else None
        if fallback is None or fallback.is_empty or fallback.geom_type != "Polygon":
            fallback = rigid.get(i, _g)
        parts[i] = (props, fallback, (note + ",rigid-clipped").strip(","))
    for i in sorted(bad):
        # a fallback plot must not sit on a sibling that received a filled piece: trim it
        props, g, note = parts[i]
        for j, (_pj, gj, _nj) in enumerate(parts):
            if j == i or gj is None or gj.is_empty or g is None:
                continue
            if g.intersection(gj).area > SIBLING_OVERLAP_M2:
                trimmed = _largest(g.difference(gj))
                if trimmed is not None and not trimmed.is_empty and trimmed.geom_type == "Polygon":
                    g = trimmed
        parts[i] = (props, g, note)
    if bad:
        notes.append("%d plot(s) rigid-clipped" % len(bad))
    return parts, notes


def report_rail_conflicts(village, report):
    p = paths.vector_dir(village) / "rail_conflicts.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["survey_a", "survey_b", "overlap_sqm", "note"])
        for a, b, area in report.get("rail_conflicts", []):
            w.writerow([a, b, area,
                        "railway land carved out of the older parcel sheet; decide in review"])
    return p
