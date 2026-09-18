"""Sheet polygons in sheet metres, and the outline geometry derived from them."""
import json
import math

import numpy as np
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union

from . import paths


def load_sheet(village, survey):
    """[(properties, Polygon)] for one sheet, in sheet metres."""
    g = json.loads(paths.sheet_path(village, survey).read_text(encoding="utf-8"))
    return [(f["properties"], shape(f["geometry"])) for f in g["features"]]


def dissolve(polys):
    """Largest single polygon of the union of a sheet's parts."""
    u = unary_union([p.buffer(0) for _props, p in polys])
    return max(u.geoms, key=lambda q: q.area) if u.geom_type == "MultiPolygon" else u


def outline(polys):
    """Exterior ring of the dissolved sheet, clockwise, every real vertex kept.

    The 2026-09-17 prototype merged vertices whose turn angle was under 1 degree. That merge moved
    vertices differently on the two sheets of a pair and broke boundary matching (spec 6.1), so
    only exact duplicates are dropped here.
    """
    ring = list(dissolve(polys).exterior.coords)[:-1]
    out = [ring[0]]
    for p in ring[1:]:
        if math.dist(p, out[-1]) > 1e-9:
            out.append(p)
    if len(out) > 1 and math.dist(out[0], out[-1]) <= 1e-9:
        out.pop()
    if not Polygon(out).exterior.is_ccw:
        return out
    return out[::-1]


def edge_lengths(v):
    return [math.dist(v[i], v[(i + 1) % len(v)]) for i in range(len(v))]


def turn_angles(v):
    """Turn angle in degrees at each vertex of a closed ring (0 = straight)."""
    n = len(v)
    out = []
    for i in range(n):
        a, b, c = v[i - 1], v[i], v[(i + 1) % n]
        t = math.degrees(math.atan2(c[1] - b[1], c[0] - b[0]) - math.atan2(b[1] - a[1], b[0] - a[0]))
        out.append(abs((t + 180) % 360 - 180))
    return out


def sample_outline(v, step=1.0):
    """Points every `step` metres along the ring, with each sample's edge bearing in radians."""
    pts, brg = [], []
    for i in range(len(v)):
        a = np.asarray(v[i], float)
        b = np.asarray(v[(i + 1) % len(v)], float)
        length = float(np.linalg.norm(b - a))
        n = max(1, int(round(length / step)))
        ang = math.atan2(b[1] - a[1], b[0] - a[0])
        for k in range(n):
            pts.append(a + (b - a) * (k + 0.5) / n)
            brg.append(ang)
    return np.array(pts), np.array(brg)
