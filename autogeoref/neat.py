"""Make a village's parcels a neat fabric: one line per seam, no overlaps, no slivers, no gaps.

    python -m autogeoref.neat <merged.geojson> <out.geojson>

Akash's ask on 2026-09-24, with his own hand-edited example at the Thirukatchur junction
(47B/54/88A/581/583): the parcels keep their shapes, and wherever two of them meet, the two
drawn boundaries become one clean line. This is the only kind of change made:

  * an overlap between two plots is split along its middle - each gives up half;
  * a gap between two plots is split along its middle - each takes half:
      - between surveys the PDF sheets print as neighbours, wedges up to 2 x GAP_PDF_M wide
        (two parcels that touch at a corner and splay apart are one seam on the PDF);
      - between any two plots, hairlines up to 2 x GAP_ANY_M wide (a double line);
  * a real lane or unplaced ground between parcels that are not neighbours is left open;
  * plots away from every seam keep their exact vertices.

Nothing is snapped wholesale and no drawing is re-fitted: the change to each plot is bounded by
half the seam width, and every change is reported.
"""
import argparse
import math

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import MultiPoint, Polygon, MultiPolygon
from shapely.strtree import STRtree

GAP_PDF_M = 3.0         # half-width of a wedge closed between PDF neighbours
GAP_ANY_M = 0.6         # half-width of a hairline closed between any two plots
MIN_PIECE_M2 = 0.05
STEP_M = 0.25           # boundary sampling for the middle-line split
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


def _points_near(line, region, reach):
    """Points every STEP_M along `line` that lie within `reach` of `region`."""
    zone = region.buffer(reach)
    part = line.intersection(zone)
    pts = []
    for seg in getattr(part, "geoms", [part]):
        if seg.is_empty or seg.length == 0:
            continue
        n = max(2, int(math.ceil(seg.length / STEP_M)) + 1)
        pts.extend((p.x, p.y) for p in (seg.interpolate(t) for t in np.linspace(0, seg.length, n)))
    return pts


def split_middle(region, side_a, side_b):
    """Split `region` along the middle line between boundary `side_a` and boundary `side_b`.

    Returns (part nearer side_a, part nearer side_b). Where the two boundaries run together they
    already are one line and say nothing about the middle, so only their separate runs count."""
    own_a = side_a.difference(side_b.buffer(0.01))
    own_b = side_b.difference(side_a.buffer(0.01))
    pa = _points_near(own_a, region, STEP_M * 2)
    pb = _points_near(own_b, region, STEP_M * 2)
    if not pa and not pb:
        return region, shapely.Polygon()
    if not pb:
        return region, shapely.Polygon()
    if not pa:
        return shapely.Polygon(), region
    # where the two boundaries meet, a point lies on both: it belongs to neither side
    from scipy.spatial import cKDTree
    pa = list(dict.fromkeys((round(x, 4), round(y, 4)) for x, y in pa))
    pb = list(dict.fromkeys((round(x, 4), round(y, 4)) for x, y in pb))
    d, _ = cKDTree(np.array(pa)).query(np.array(pb))
    pb = [p for p, dd in zip(pb, d) if dd > 1e-3]
    if not pb:
        return region, shapely.Polygon()
    pts = MultiPoint(pa + pb)
    cells = shapely.voronoi_polygons(pts, extend_to=region.buffer(10), ordered=True)
    cells = list(cells.geoms)
    near_a = shapely.union_all(cells[:len(pa)])
    return _inter(region, near_a), _diff(region, near_a)


def split_nearest(region, sides):
    """Split `region` among several plots: each part goes to the plot whose boundary is nearest.

    `sides` is {key: boundary}. Points lying on more than one boundary are dropped."""
    from collections import Counter
    pts, keys = [], []
    for k, line in sides.items():
        for x, y in _points_near(line, region, STEP_M * 2):
            pts.append((round(x, 4), round(y, 4)))
            keys.append(k)
    count = Counter(pts)
    keep = [(p, k) for p, k in zip(pts, keys) if count[p] == 1]
    keep = list(dict(keep).items())
    if len({k for _p, k in keep}) < 2:
        return {keep[0][1]: region} if keep else {}
    cells = list(shapely.voronoi_polygons(MultiPoint([p for p, _k in keep]), extend_to=region.buffer(10),
                                          ordered=True).geoms)
    out = {}
    for k in {k for _p, k in keep}:
        part = _inter(region, shapely.union_all([c for c, (_p, kk) in zip(cells, keep) if kk == k]))
        if not part.is_empty:
            out[k] = part
    return out


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


def _inter(a, b):
    return _op(shapely.intersection, a, b)


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


def clean(gdf, printed, survey_col="survey_no"):
    """(clean GeoDataFrame, report DataFrame). `printed`: {survey: set of PDF-neighbour surveys}."""
    original = list(gdf.geometry)
    geoms = [_clean_geom(g) if g is not None else shapely.Polygon() for g in original]
    surveys = [str(s) for s in gdf[survey_col]]
    before = list(geoms)
    touched = set()

    def neighbours_pdf(i, j):
        a, b = surveys[i], surveys[j]
        return a != b and (b in printed.get(a, set()) or a in printed.get(b, set()))

    # 1. overlaps: split along their middle
    for _round in range(3):
        tree = STRtree(geoms)
        found = 0
        for i in range(len(geoms)):
            for j in tree.query(geoms[i]):
                if j <= i:
                    continue
                inter = _inter(geoms[i], geoms[j])
                for piece in _area_pieces(inter, 1e-3):
                    # the half next to i's own edge (inside j) goes to j, and vice versa
                    to_j, to_i = split_middle(piece, geoms[i].boundary, geoms[j].boundary)
                    geoms[i] = _clean_geom(_diff(geoms[i], to_j))
                    geoms[j] = _clean_geom(_diff(geoms[j], to_i))
                    touched.update((i, j))
                    found += 1
        if not found:
            break

    # 2. gaps: wedges between PDF neighbours, hairlines between any plots
    tree = STRtree(geoms)
    for i in range(len(geoms)):
        for j in tree.query(geoms[i].buffer(2 * GAP_PDF_M)):
            if j <= i:
                continue
            r = GAP_PDF_M if neighbours_pdf(i, j) else GAP_ANY_M
            a, b = geoms[i], geoms[j]
            if a.is_empty or b.is_empty or a.distance(b) > 2 * r:
                continue
            u = _union(a, b)
            closed = u.buffer(r, join_style=2, mitre_limit=2).buffer(-r, join_style=2, mitre_limit=2)
            wedge = _inter(_inter(_diff(closed, u), a.buffer(r + 0.05)), b.buffer(r + 0.05))
            if wedge.is_empty:
                continue
            others = [k for k in tree.query(wedge) if k not in (i, j)]
            if others:
                wedge = _diff(wedge, shapely.union_all([geoms[k] for k in others], grid_size=GRID_M))
            for piece in _area_pieces(wedge):
                if piece.distance(a) > 0.01 or piece.distance(b) > 0.01:
                    continue
                to_a, to_b = split_middle(piece, a.boundary, b.boundary)
                geoms[i] = _clean_geom(_union(geoms[i], to_a))
                geoms[j] = _clean_geom(_union(geoms[j], to_b))
                a, b = geoms[i], geoms[j]
                touched.update((i, j))

    # 2b. wedges between PDF-neighbour surveys where three or more plots meet
    by_survey = {}
    for i, s_ in enumerate(surveys):
        by_survey.setdefault(s_, []).append(i)
    done = set()
    for sa, ia in by_survey.items():
        for sb in sorted(printed.get(sa, set())):
            if sb not in by_survey or (sb, sa) in done:
                continue
            done.add((sa, sb))
            A = shapely.union_all([geoms[i] for i in ia], grid_size=GRID_M)
            B = shapely.union_all([geoms[i] for i in by_survey[sb]], grid_size=GRID_M)
            if A.is_empty or B.is_empty or A.distance(B) > 2 * GAP_PDF_M:
                continue
            u = _union(A, B)
            closed = u.buffer(GAP_PDF_M, join_style=2, mitre_limit=2).buffer(-GAP_PDF_M, join_style=2, mitre_limit=2)
            wedge = _inter(_inter(_diff(closed, u), A.buffer(GAP_PDF_M + 0.05)), B.buffer(GAP_PDF_M + 0.05))
            if wedge.is_empty:
                continue
            tree = STRtree(geoms)
            near = list(tree.query(wedge.buffer(0.05)))
            others = [k for k in near if surveys[k] not in (sa, sb)]
            if others:
                wedge = _diff(wedge, shapely.union_all([geoms[k] for k in others], grid_size=GRID_M))
            for piece in _area_pieces(wedge):
                if piece.distance(A) > 0.01 or piece.distance(B) > 0.01:
                    continue
                cand = {k: geoms[k].boundary for k in near if surveys[k] in (sa, sb) and geoms[k].distance(piece) < 0.05}
                for k, part in split_nearest(piece, cand).items():
                    geoms[k] = _clean_geom(_union(geoms[k], part))
                    touched.add(k)

    # 3. make every shared line exactly one line: snap touched plots onto their neighbours
    tree = STRtree(geoms)
    for i in sorted(touched):
        near = [k for k in tree.query(geoms[i].buffer(0.05)) if k != i]
        if near:
            ref = shapely.union_all([geoms[k].boundary for k in near])
            geoms[i] = _clean_geom(shapely.snap(geoms[i], ref, 0.02))
    # any overlap left (a snap, or a hairline under the overlap threshold) is removed from the
    # later plot; every plot changed here is reported and written
    for _round in range(2):
        tree = STRtree(geoms)
        for i in range(len(geoms)):
            for j in tree.query(geoms[i]):
                if j > i and _inter(geoms[i], geoms[j]).area > 1e-4:
                    geoms[j] = _clean_geom(_diff(geoms[j], geoms[i]))
                    touched.add(j)

    out = gdf.copy()
    out["geometry"] = [g if i in touched else original[i] for i, g in enumerate(geoms)]
    out = out.set_geometry("geometry")
    rep = pd.DataFrame({
        "row": range(len(out)), "survey_no": surveys,
        "area_before_m2": [round(g.area, 2) for g in before],
        "area_after_m2": [round(g.area, 2) for g in out.geometry],
        "moved_max_m": [round(before[i].hausdorff_distance(out.geometry.iloc[i]), 2) if i in touched else 0.0
                        for i in range(len(out))],
        "edited": [i in touched for i in range(len(out))]})
    rep["area_change_m2"] = (rep.area_after_m2 - rep.area_before_m2).round(2)
    return out, rep


def main(argv=None):
    import geopandas as gpd
    from . import neighbours
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
    junk = out.geometry.is_empty | out.geometry.isna() | ((out.geometry.area < 0.05) & out.get("plot_no", pd.Series(index=out.index)).isna())
    rep["dropped_as_junk"] = junk.values
    out = out[~junk].copy()
    # his mandatory fields first: kide / survey_no / subdiv_no / area_acre, then area_sqm
    from . import fmb_finish
    keys = [fmb_finish.survey_keys(s_, pl) for s_, pl in zip(out.survey_no.astype(str), out.get("plot_no", [None] * len(out)))]
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
