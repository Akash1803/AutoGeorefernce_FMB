"""Make a village's parcels a neat fabric: one line per seam, no overlaps, no slivers, no gaps.

    python -m autogeoref.neat <merged.geojson> <out.geojson>

Akash's ask on 2026-09-24, with his own hand-edited example at the Thirukatchur junction
(47B/54/88A/581/583): parcels keep their shapes, and where two of them meet, the gap or overlap
is closed onto the nearby parcel's existing line - the way it is done by hand with snapping:

  * parcels are settled largest first; a later parcel's vertex near a settled neighbour moves
    onto that neighbour's existing corner (corners meet corners) or edge - up to GAP_PDF_M between
    surveys the PDF sheets print as neighbours, up to GAP_ANY_M between any two plots;
  * every plot sharing a moved vertex moves with it, so a survey's own subdivision lines stay joined;
  * a neighbour's corner lying on the moved edge is added, and what still overlaps a settled
    neighbour is clipped along the neighbour's existing line;
  * a real lane or unplaced ground between parcels that are not neighbours is left open;
  * plots away from every seam keep their exact vertices.

No new boundary is drawn: every line in the result was already in the data. (A first version split
gaps along their middle; Akash rejected the curved lines it drew inside the gaps.)
"""
import argparse
import math

import numpy as np
import pandas as pd
import shapely
import shapely.ops
from scipy.spatial import cKDTree
from shapely.geometry import MultiPolygon, Polygon
from shapely.strtree import STRtree

GAP_PDF_M = 5.0         # how far a vertex may move onto a PDF neighbour (never through a parcel)
GAP_ANY_M = 0.6         # how far a vertex may move onto any plot (a double line)
WEDGE_M = 1.5           # a wedge left between PDF neighbours is closed up to 2 x this width
WEDGE_MAX_M2 = 40.0     # ... and only up to this size; anything bigger is open ground
MIN_PIECE_M2 = 0.05
GRID_M = 0.001


def _polys(g):
    if g is None or g.is_empty:
        return []
    if isinstance(g, Polygon):
        return [g]
    if isinstance(g, MultiPolygon):
        return list(g.geoms)
    return [p for x in getattr(g, "geoms", []) for p in _polys(x)]


def _area_pieces(g, min_area=MIN_PIECE_M2):
    return [p for p in _polys(g) if p.area > min_area]


def _op(fn, a, b):
    """Overlay on the GRID_M grid (robust); on a GEOS topology error, repair and retry."""
    try:
        return fn(a, b, grid_size=GRID_M)
    except shapely.errors.GEOSException:
        return fn(shapely.make_valid(a).buffer(0), shapely.make_valid(b).buffer(0), grid_size=GRID_M)


def _diff(a, b):
    return _op(shapely.difference, a, b)


def _union(a, b):
    return _op(shapely.union, a, b)


def _clean_geom(g):
    g = shapely.make_valid(g)
    try:
        g = shapely.make_valid(shapely.set_precision(g, GRID_M))
    except shapely.errors.GEOSException:
        g = shapely.make_valid(shapely.set_precision(g.buffer(0), GRID_M))
    parts = _area_pieces(g, 1e-4)
    if not parts:
        return shapely.Polygon()
    return parts[0] if len(parts) == 1 else MultiPolygon(parts)


def _key(x, y):
    return (round(x, 3), round(y, 3))


def _move_nodes(g, moves):
    """Apply {node: new xy} to every vertex of g that sits on a moved node."""
    if not moves or g.is_empty:
        return g

    def f(coords):
        out = coords.copy()
        for n, (x, y) in enumerate(coords):
            t = moves.get(_key(x, y))
            if t is not None:
                out[n] = t
        return out
    return shapely.transform(g, f)


def _snap_target(pt, lines, nodes, tol):
    """Where a vertex closes onto a neighbour: its nearest corner when that is about as near as
    the nearest edge (corners meet corners), else the nearest point on the edge; None beyond tol."""
    d_edge = lines.distance(pt)
    if d_edge > tol:
        return None
    if nodes is not None and len(nodes):
        d_node, k = cKDTree(nodes).query([pt.x, pt.y])
        if d_node <= tol and d_node <= d_edge + 0.2:           # never slide a vertex along a line
            return (float(nodes[k][0]), float(nodes[k][1]))
    q = shapely.ops.nearest_points(lines, pt)[0]
    return (q.x, q.y)


def _internal_direction(node, i, geoms, surveys, tree):
    """If `node` is where one of this survey's internal lines meets its outline, the unit
    direction of that line arriving at the node (from inside); else None."""
    siblings = [k for k in tree.query(shapely.Point(node).buffer(0.002))
                if surveys[k] == surveys[i] and not geoms[k].is_empty]
    if len(siblings) < 2:
        return None
    neighbours_of = {}
    for k in siblings:
        for part in _polys(geoms[k]):
            c = np.asarray(part.exterior.coords)[:-1]
            keys = [_key(x, y) for x, y in c]
            for n, kk in enumerate(keys):
                if kk == node:
                    for m in ((n - 1) % len(c), (n + 1) % len(c)):
                        neighbours_of.setdefault(keys[m], set()).add(k)
    shared = [kk for kk, ks in neighbours_of.items() if len(ks) >= 2]  # the internal edge's other end
    if len(shared) != 1:
        return "fixed"                     # several internal lines meet here: it never moves sideways
    d = np.array(node) - np.array(shared[0])
    n = np.hypot(*d)
    return d / n if n > 1e-6 else None


def _crosses_a_parcel(v, t, i, geoms, tree):
    """True when moving vertex v to t would pass through a parcel other than the one it closes
    onto (e.g. a big plot snapping across a narrow strip to the plot beyond it)."""
    path = shapely.LineString([v, t])
    for k in tree.query(path):
        if k == i or geoms[k].is_empty:
            continue
        inside = geoms[k].buffer(-0.02)
        if not inside.is_empty and path.intersection(inside).length > 0.05:
            return True
    return False


def _acceptable(new, old):
    if new.is_empty or not new.is_valid or len(_polys(new)) > max(1, len(_polys(old))):
        return False
    return new.area >= 0.85 * old.area


def clean(gdf, printed, survey_col="survey_no", fixed_surveys=()):
    """(clean GeoDataFrame, report DataFrame). `printed`: {survey: set of PDF-neighbour surveys}.

    `fixed_surveys` are the base: never changed in any way (Akash's hand-placed rail strip 582,
    parcels he has approved); everything else closes onto them."""
    original = list(gdf.geometry)
    geoms = [_clean_geom(g) if g is not None else shapely.Polygon() for g in original]
    surveys = [str(s) for s in gdf[survey_col]]
    before = list(geoms)
    touched = set()

    def pdf_pair(i, j):
        a, b = surveys[i], surveys[j]
        return a != b and (b in printed.get(a, set()) or a in printed.get(b, set()))

    frozen = {i for i, s_ in enumerate(surveys) if s_ in {str(x) for x in fixed_surveys}}
    settled = set(frozen)
    for i in sorted(range(len(geoms)), key=lambda i: -geoms[i].area):
        if i in frozen:
            continue
        g = geoms[i]
        if g.is_empty:
            settled.add(i)
            continue
        tree = STRtree(geoms)
        near = [k for k in tree.query(g.buffer(2 * GAP_PDF_M))
                if k != i and k in settled and not geoms[k].is_empty]
        if not near:
            settled.add(i)
            continue
        fixed = {_key(x, y) for k in near for x, y in shapely.get_coordinates(geoms[k])}
        near_union = shapely.union_all([geoms[k] for k in near], grid_size=GRID_M)
        ov0 = _op(shapely.intersection, g, near_union).area          # overlap it really had
        allowed = ov0 + max(1.0, 0.02 * g.area)
        pdf_refs = [k for k in near if pdf_pair(i, k)]
        done = None
        for scale in (1.0, 0.5, 0.25, 0.0):
            moves = {}
            for refs, tol in ((pdf_refs, GAP_PDF_M * scale), (near, GAP_ANY_M * min(1.0, scale * 4))):
                if not refs or tol <= 0:
                    continue
                lines = shapely.union_all([geoms[k].boundary for k in refs])
                nodes = np.unique(np.vstack([shapely.get_coordinates(geoms[k]) for k in refs]), axis=0)
                for x, y in shapely.get_coordinates(g):
                    kk = _key(x, y)
                    if kk in fixed or kk in moves:
                        continue
                    t = _snap_target(shapely.Point(x, y), lines, nodes, tol)
                    if t is None or math.dist(t, (x, y)) <= 1e-4:
                        continue
                    owners = {surveys[k] for k in tree.query(shapely.Point(x, y).buffer(0.002))
                              if not geoms[k].is_empty and geoms[k].boundary.distance(shapely.Point(x, y)) < 0.002}
                    if len(owners) >= 2 and math.dist(t, (x, y)) > 0.05:
                        continue                  # a corner two surveys already share never moves
                    d = _internal_direction(kk, i, geoms, surveys, tree)
                    if isinstance(d, str):
                        if math.dist(t, (x, y)) > 0.05:
                            continue
                    elif d is not None and math.dist(t, (x, y)) > 0.05:
                        # an internal subdivision line is extended (or shortened) along itself
                        ray = shapely.LineString([(x - d[0] * tol, y - d[1] * tol), (x + d[0] * tol, y + d[1] * tol)])
                        hit = ray.intersection(lines)
                        pts = [q for q in getattr(hit, "geoms", [hit]) if q.geom_type == "Point"]
                        if not pts:
                            continue
                        q = min(pts, key=lambda q: q.distance(shapely.Point(x, y)))
                        t = (q.x, q.y)
                    if _crosses_a_parcel((x, y), t, i, geoms, tree):
                        continue                                  # never snap across another parcel
                    moves[kk] = t
            trial = {j: _clean_geom(_move_nodes(geoms[j], moves)) for j in range(len(geoms))
                     if moves and j not in settled and not geoms[j].is_empty
                     and any(_key(x, y) in moves for x, y in shapely.get_coordinates(geoms[j]))}
            if not all(_acceptable(trial[j], geoms[j]) for j in trial):
                continue
            gi = trial.get(i, g)
            # closing a gap must not push this plot over a third plot that is still to be settled
            pending = [k for k in tree.query(gi) if k != i and k not in settled and not geoms[k].is_empty]
            if pending:
                pend = shapely.union_all([geoms[k] for k in pending], grid_size=GRID_M)
                if _op(shapely.intersection, gi, pend).area - _op(shapely.intersection, g, pend).area > 0.5:
                    continue
            # a settled neighbour's corner lying on the moved edge joins it; then clip along its line
            ref = shapely.union_all([geoms[k].boundary for k in near])
            g3 = _clean_geom(_diff(_clean_geom(shapely.snap(gi, ref, 0.05)), near_union))
            if g3.is_empty or g.area - g3.area > allowed:
                continue                                              # it would lose real ground
            parts = _polys(g3)
            if len(parts) > max(1, len(_polys(before[i]))):
                body = max(parts, key=lambda q: q.area)
                if g3.area - body.area > 1.0:
                    continue                                          # a real piece would be cut off
                g3 = body
            done = (trial, g3)
            break
        if done is not None:
            trial, g3 = done
            for j, gj in trial.items():
                if j != i and not gj.equals(geoms[j]):
                    geoms[j] = gj
                    touched.add(j)
            if not g3.equals(geoms[i]):
                geoms[i] = g3
                touched.add(i)
        settled.add(i)

    # where a moved vertex now lies on a neighbour's edge, the neighbour gets that vertex too, so
    # the seam is one shared line vertex for vertex (its shape does not change)
    tree = STRtree(geoms)
    for i in sorted(touched):
        for k in tree.query(geoms[i].buffer(0.01)):
            if k == i or k in frozen or geoms[k].is_empty:
                continue
            snapped = shapely.snap(geoms[k], geoms[i], 0.002)
            if shapely.get_num_coordinates(snapped) != shapely.get_num_coordinates(geoms[k]) \
                    and snapped.is_valid and abs(snapped.area - geoms[k].area) < 0.01:
                geoms[k] = snapped
                touched.add(k)

    # hairline slivers still open between plots: each joins the plot it shares the longest edge with
    u = shapely.union_all(geoms, grid_size=GRID_M)
    closed = u.buffer(GAP_ANY_M, join_style=2, mitre_limit=2).buffer(-GAP_ANY_M, join_style=2, mitre_limit=2)
    tree = STRtree(geoms)
    for piece in _area_pieces(_diff(closed, u), 0.01):
        cand = [k for k in tree.query(piece.buffer(0.01)) if not geoms[k].is_empty]
        if len(cand) < 2:
            continue
        cand = [k for k in cand if geoms[k].area >= 4 * piece.area and k not in frozen]
        if not cand:
            continue
        k = max(cand, key=lambda k: piece.boundary.intersection(geoms[k].buffer(0.01)).length)
        geoms[k] = _clean_geom(_union(geoms[k], piece))
        touched.add(k)

    # a wedge still open between PDF neighbours (too far to snap safely) joins the neighbouring
    # plot it shares the longest edge with, so the seam closes on that plot's existing line
    by_survey = {}
    for i, s_ in enumerate(surveys):
        by_survey.setdefault(s_, []).append(i)
    seen = set()
    for sa, ia in by_survey.items():
        for sb in sorted(printed.get(sa, set())):
            if sb not in by_survey or (sb, sa) in seen:
                continue
            seen.add((sa, sb))
            A = shapely.union_all([geoms[i] for i in ia], grid_size=GRID_M)
            B = shapely.union_all([geoms[i] for i in by_survey[sb]], grid_size=GRID_M)
            if A.is_empty or B.is_empty or A.distance(B) > 0.05:
                continue
            u = _union(A, B)
            closed = u.buffer(WEDGE_M, join_style=2, mitre_limit=2).buffer(-WEDGE_M, join_style=2, mitre_limit=2)
            wedge = _diff(closed, u)
            if wedge.is_empty:
                continue
            tree = STRtree(geoms)
            others = [k for k in tree.query(wedge) if surveys[k] not in (sa, sb) and not geoms[k].is_empty]
            if others:
                wedge = _diff(wedge, shapely.union_all([geoms[k] for k in others], grid_size=GRID_M))
            for piece in _area_pieces(wedge, 0.3):
                if piece.distance(A) > 0.01 or piece.distance(B) > 0.01 or piece.area > WEDGE_MAX_M2:
                    continue                                 # open ground, not a seam: leave it
                cand = [k for k in tree.query(piece.buffer(0.01)) if surveys[k] in (sa, sb) and k not in frozen
                        and not geoms[k].is_empty and geoms[k].area >= 4 * piece.area]
                if not cand:
                    continue
                k = max(cand, key=lambda k: piece.boundary.intersection(geoms[k].buffer(0.01)).length)
                geoms[k] = _clean_geom(_union(geoms[k], piece))
                touched.add(k)

    # the base wins: whatever still overlaps a fixed parcel is trimmed along the fixed parcel's
    # line; any overlap left between two other plots goes to the larger one
    if frozen:
        base = shapely.union_all([geoms[k] for k in frozen if not geoms[k].is_empty], grid_size=GRID_M)
        for i in range(len(geoms)):
            if i in frozen or geoms[i].is_empty or _op(shapely.intersection, geoms[i], base).area < 0.05:
                continue
            cut = _clean_geom(_diff(geoms[i], base))
            if not cut.is_empty:
                geoms[i] = cut
                touched.add(i)
    tree = STRtree(geoms)
    for i in range(len(geoms)):
        if i in frozen or geoms[i].is_empty:
            continue
        for j in tree.query(geoms[i]):
            if j == i or j in frozen or geoms[j].is_empty:
                continue
            if _op(shapely.intersection, geoms[i], geoms[j]).area > 0.05:
                small, big = (i, j) if geoms[i].area < geoms[j].area else (j, i)
                geoms[small] = _clean_geom(_diff(geoms[small], geoms[big]))
                touched.add(small)

    out = gdf.copy()
    out["geometry"] = [g if i in touched else original[i] for i, g in enumerate(geoms)]
    out = out.set_geometry("geometry")
    rep = pd.DataFrame({
        "row": range(len(out)), "survey_no": surveys,
        "area_before_m2": [round(g.area, 2) for g in before],
        "area_after_m2": [round(g.area, 2) for g in out.geometry],
        "moved_max_m": [round(before[i].hausdorff_distance(out.geometry.iloc[i]), 2)
                        if i in touched and not out.geometry.iloc[i].is_empty and not before[i].is_empty else 0.0
                        for i in range(len(out))],
        "edited": [i in touched for i in range(len(out))]})
    rep["area_change_m2"] = (rep.area_after_m2 - rep.area_before_m2).round(2)
    return out, rep


def main(argv=None):
    import geopandas as gpd
    from . import fmb_finish, neighbours
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("src")
    ap.add_argument("out")
    a = ap.parse_args(argv)
    g = gpd.read_file(a.src)
    crs = g.crs
    g = g.to_crs(32644)
    # a merged layer can carry rows without village fields (approved _auto layers): take the
    # village from the file each row came from, and its name/district from its village-mates
    if "path" in g.columns and "village_code" in g.columns:
        from_path = g["path"].astype(str).str.replace("\\", "/").str.extract(r"FMB_Vector/(\d+_\d+_\d+)/")[0]
        g["village_code"] = g["village_code"].fillna(from_path)
        for col in ("village", "district", "district_code", "taluk", "taluk_code"):
            if col in g.columns:
                g[col] = g[col].fillna(g.groupby("village_code")[col].transform("first"))
    printed = {}
    for v in sorted(set(g.get("village_code", pd.Series(dtype=str)).dropna().astype(str))):
        for s in set(g.loc[g.village_code == v, "survey_no"].astype(str)):
            printed.setdefault(s, set()).update(neighbours.printed(v, s))
    out, rep = clean(g, printed)
    # junk rows (no geometry, or a hairline with no plot number) do not belong in a ready file
    junk = out.geometry.is_empty | out.geometry.isna() | (
        (out.geometry.area < 0.05) & out.get("plot_no", pd.Series(index=out.index)).isna())
    rep["dropped_as_junk"] = junk.values
    out = out[~junk].copy()
    # his mandatory fields first: kide / survey_no / subdiv_no / area_acre, then area_sqm
    keys = [fmb_finish.survey_keys(s_, pl) for s_, pl in
            zip(out.survey_no.astype(str), out.get("plot_no", [None] * len(out)))]
    out["kide"] = [k for _sub, k in keys]
    out["subdiv_no"] = [sub for sub, _k in keys]
    out["area_sqm"] = out.geometry.area.round(2)
    out["area_acre"] = (out.geometry.area / 4046.8564224).round(4)
    lead = ["kide", "survey_no", "subdiv_no", "area_acre", "area_sqm"]
    out = out[lead + [c for c in out.columns if c not in lead]]
    out.to_crs(crs).to_file(a.out, driver="GeoJSON")
    rep.to_csv(str(a.out).rsplit(".", 1)[0] + "_changes.csv", index=False)
    print("edited %d of %d plots; largest move %.2f m; dropped %d junk rows"
          % (rep.edited.sum(), len(rep), rep.moved_max_m.max(), int(junk.sum())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
