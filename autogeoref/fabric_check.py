"""The one topology and shape check, run on the parcels exactly as Akash's QGIS shows them.

    python -m autogeoref.fabric_check 35_04_074 [--csv <folder>]

Written on 2026-09-24 after the audit found that every "clean" report of 23-24 Sep had measured
hidden copies with minimum distance. This check:

  * reads each survey through visible.shown()/read_hand() - the file and layer QGIS draws;
  * overlaps: area shared by two surveys (> 0.5 m2);
  * gaps: only between surveys the PDF sheets print as neighbours, measured as the open wedge
    along their shared edge - two parcels that touch at one corner and splay apart are a gap,
    even though their minimum distance is 0.00 m;
  * shape: the PDF drawing laid rigidly (no scale, no mirror) on the parcel; IoU 1.00 = identical.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
from shapely import affinity
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely.strtree import STRtree
from scipy.spatial import cKDTree

from . import neighbours, paths, sheets, visible

OVERLAP_MIN_M2 = 0.5
GAP_MIN_M2 = 1.0
WEDGE_REACH_M = 3.0
FAR_M = 10.0


def _polys(g):
    if g is None or g.is_empty:
        return []
    if isinstance(g, Polygon):
        return [g]
    if isinstance(g, MultiPolygon):
        return list(g.geoms)
    return [p for x in getattr(g, "geoms", []) for p in _polys(x)]


def overlaps(bodies, min_m2=OVERLAP_MIN_M2):
    """[(a, b, m2)] for every pair of surveys that claim the same ground."""
    names = list(bodies)
    geoms = [bodies[n] for n in names]
    tree = STRtree(geoms)
    out = []
    for i, g in enumerate(geoms):
        for j in tree.query(g):
            if j > i:
                a = g.intersection(geoms[j]).area
                if a > min_m2:
                    out.append((names[i], names[j], round(a, 2)))
    return out


def seam_gap(a, b, reach=WEDGE_REACH_M, others=None):
    """Open area between two parcels along their shared edge (m2). Ground covered by any other
    parcel (`others`, a geometry) is not open."""
    u = unary_union([a, b])
    wedge = (u.buffer(reach, join_style=2).buffer(-reach, join_style=2).difference(u)
             .intersection(a.buffer(reach + 0.2)).intersection(b.buffer(reach + 0.2)))
    if others is not None and not others.is_empty and not wedge.is_empty:
        wedge = wedge.difference(others)
    return sum(p.area for p in _polys(wedge) if p.area > 0.3)


def gaps(bodies, printed, min_m2=GAP_MIN_M2):
    """[(a, b, m2, note)] open gaps between surveys the PDF sheets name as neighbours."""
    out, done = [], set()
    tree_names = list(bodies)
    tree = STRtree([bodies[n] for n in tree_names])
    for a in sorted(bodies):
        for b in sorted(printed.get(a, set())):
            if b not in bodies or (b, a) in done or (a, b) in done:
                continue
            done.add((a, b))
            d = bodies[a].distance(bodies[b])
            if d > FAR_M:
                out.append((a, b, None, "PDF neighbours placed %.1f m apart" % d))
                continue
            zone = bodies[a].union(bodies[b]).buffer(WEDGE_REACH_M + 0.5)
            near = [bodies[tree_names[k]] for k in tree.query(zone) if tree_names[k] not in (a, b)]
            m2 = seam_gap(bodies[a], bodies[b], others=unary_union(near) if near else None)
            if m2 > min_m2:
                out.append((a, b, round(m2, 1), "touching at a point, open along the edge" if d < 0.05 else "apart %.2f m" % d))
    return out


def _sample(geom, step=0.75):
    pts = []
    for p in _polys(geom):
        for ring in [p.exterior, *p.interiors]:
            n = max(3, int(ring.length / step))
            pts.extend((q.x, q.y) for q in (ring.interpolate(t) for t in np.linspace(0, ring.length, n)))
    return np.array(pts)


def _kabsch(a, b):
    ca, cb = a.mean(0), b.mean(0)
    u, _s, vt = np.linalg.svd((a - ca).T @ (b - cb))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, d]) @ u.T
    return r, cb - r @ ca


def _apply(g, r, t):
    return affinity.affine_transform(g, [r[0, 0], r[0, 1], r[1, 0], r[1, 1], t[0], t[1]])


def _iou(a, b):
    u = a.union(b).area
    return a.intersection(b).area / u if u else 0.0


def shape_iou(drawing, placed):
    """(IoU, (R, t)): the PDF drawing laid rigidly on the placed parcel.

    Point-to-point ICP from 36 starting headings finds the basin; it slides slowly along straight
    edges, so the best three are finished by maximising the overlap itself (Nelder-Mead over
    shift and heading). No scale, no mirror.
    """
    from scipy.optimize import minimize
    src, dst = _sample(drawing), _sample(placed)
    if len(src) < 3 or len(dst) < 3:
        return float("nan"), None
    tree = cKDTree(dst)
    cs, cd = np.array(drawing.centroid.coords[0]), np.array(placed.centroid.coords[0])
    starts = []
    for deg in range(0, 360, 10):
        a = math.radians(deg)
        r = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
        t = cd - r @ cs
        for _ in range(40):
            d, i = tree.query(src @ r.T + t)
            keep = d < max(3.0, np.median(d) * 2.5)
            if keep.sum() < 3:
                break
            r, t = _kabsch(src[keep], dst[i[keep]])
        starts.append((_iou(_apply(drawing, r, t), placed), math.atan2(r[1, 0], r[0, 0]), t))
    starts.sort(key=lambda x: -x[0])
    centred = affinity.translate(drawing, -cs[0], -cs[1])

    def pose(x):
        return affinity.translate(affinity.rotate(centred, x[2], origin=(0, 0), use_radians=True), x[0], x[1])

    best = (-1.0, None)
    for _iou0, th, t in starts[:3]:
        c0 = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]]) @ cs + t
        res = minimize(lambda x: -_iou(pose(x), placed), [c0[0], c0[1], th], method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-5, "maxiter": 400,
                                "initial_simplex": [[c0[0], c0[1], th], [c0[0] + 1, c0[1], th],
                                                    [c0[0], c0[1] + 1, th], [c0[0], c0[1], th + 0.02]]})
        if -res.fun > best[0]:
            x = res.x
            r = np.array([[math.cos(x[2]), -math.sin(x[2])], [math.sin(x[2]), math.cos(x[2])]])
            best = (-res.fun, (r, np.array([x[0], x[1]]) - r @ cs))
    return round(best[0], 3), best[1]


def check_village(village):
    """{'overlaps': [...], 'gaps': [...], 'shapes': [...]} for the parcels as QGIS shows them."""
    bodies, shapes = {}, []
    for s, (f, _layer) in sorted(visible.shown(village).items(), key=lambda kv: paths.survey_sort_key(kv[0])):
        g = visible.read_hand(f)
        if g.crs is not None and g.crs.to_epsg() != 32644:
            g = g.to_crs(32644)
        parts = [x.buffer(0) for x in g.geometry if x is not None and not x.is_empty]
        body = unary_union(parts)
        bodies[s] = body
        row = {"survey": s, "file": f.name, "plots": len(parts), "area_m2": round(body.area, 1)}
        if paths.sheet_path(village, s).exists():
            sheet = sheets.load_sheet(village, s)
            drawing = unary_union([p.buffer(0) for _pr, p in sheet])
            row["pdf_plots"] = len(sheet)
            row["area_vs_pdf"] = round(body.area / drawing.area, 3) if drawing.area else None
            row["iou_vs_pdf"] = shape_iou(drawing, body)[0]
        shapes.append(row)
    printed = {s: {n for n in neighbours.printed(village, s)} for s in bodies}
    for s in list(printed):
        for n in printed[s]:
            printed.setdefault(n, set()).add(s)
    return {"overlaps": overlaps(bodies), "gaps": gaps(bodies, printed), "shapes": shapes}


def write_csv(village, result, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    out = folder / ("fabric_check_%s.csv" % village)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "survey", "other", "value", "note"])
        for a, b, m2 in result["overlaps"]:
            w.writerow(["overlap", a, b, m2, "m2 claimed by both"])
        for a, b, m2, note in result["gaps"]:
            w.writerow(["gap", a, b, m2 if m2 is not None else "", note])
        for r in result["shapes"]:
            w.writerow(["shape", r["survey"], "", r.get("iou_vs_pdf", ""),
                        "IoU vs PDF; area %s of PDF; %s of %s plots" % (r.get("area_vs_pdf"), r["plots"], r.get("pdf_plots"))])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("village")
    ap.add_argument("--csv", help="folder for fabric_check_<village>.csv")
    a = ap.parse_args(argv)
    res = check_village(a.village)
    print("%s: %d overlaps (%.0f m2), %d gaps between PDF neighbours, %d surveys"
          % (a.village, len(res["overlaps"]), sum(x[2] for x in res["overlaps"]),
             len(res["gaps"]), len(res["shapes"])))
    if a.csv:
        print("  written:", write_csv(a.village, res, a.csv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
