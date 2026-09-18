"""Straight edges visible in the satellite raster: walls, roads, bunds, the railway fence.

Measured on the 2026-09-18 spike: 53 % of the team's outline length lies within 1.5 m of a
same-direction segment, but 92 % of that is parallel to the corridor. Imagery therefore fixes
rotation and cross-track offset; along-track position has to come from neighbours.
"""
import math
from collections import namedtuple

import cv2
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy.spatial import cKDTree

Segment = namedtuple("Segment", "p0 p1 bearing length")


def compass(bearing_rad):
    """Segment bearings are maths angles anticlockwise from +x; this gives the compass bearing.

    The corridor here runs at compass 25.8 degrees, which is a maths angle of 64.2.
    """
    return (90.0 - math.degrees(bearing_rad)) % 180.0


def detect(raster_path, bounds=None, min_len_m=3.0):
    """Line segments in EPSG:32644 with bearings in radians 0..pi."""
    with rasterio.open(str(raster_path)) as src:
        if bounds is not None:
            win = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
            arr = src.read(window=win)
            tr = src.window_transform(win)
        else:
            arr = src.read()
            tr = src.transform
        res = src.res[0]
    if arr.size == 0:
        return []
    gray = cv2.cvtColor(np.ascontiguousarray(np.moveaxis(arr, 0, -1)), cv2.COLOR_RGB2GRAY)
    gray = cv2.bilateralFilter(gray, 7, 40, 7)
    lines = cv2.createLineSegmentDetector(0).detect(gray)[0]
    out = []
    for line in (lines if lines is not None else []):
        x0, y0, x1, y1 = [float(v) for v in np.ravel(line)[:4]]
        if math.hypot(x1 - x0, y1 - y0) * res < min_len_m:
            continue
        a = tr * (x0, y0)
        b = tr * (x1, y1)
        out.append(Segment((a[0], a[1]), (b[0], b[1]),
                           math.atan2(b[1] - a[1], b[0] - a[0]) % math.pi, math.dist(a, b)))
    return out


def directions(segments, tol_deg=10.0):
    """[(bearing, total_length)] after grouping segments whose bearings agree within tol_deg."""
    tol = math.radians(tol_deg)
    groups = []
    for s in sorted(segments, key=lambda s: -s.length):
        for i, (b, L) in enumerate(groups):
            if min(abs(s.bearing - b), math.pi - abs(s.bearing - b)) <= tol:
                groups[i] = (b, L + s.length)
                break
        else:
            groups.append((s.bearing, s.length))
    return sorted(groups, key=lambda g: -g[1])


class SegmentIndex:
    """KD-tree over points sampled along every segment, so a pose is scored in one vectorised query."""

    def __init__(self, segments, step=0.5):
        pts, brg = [], []
        for s in segments:
            a = np.array(s.p0, float)
            b = np.array(s.p1, float)
            n = max(2, int(s.length / step))
            for k in range(n + 1):
                pts.append(a + (b - a) * k / n)
                brg.append(s.bearing)
        self.points = np.array(pts) if pts else np.zeros((0, 2))
        self.bearings = np.array(brg) if brg else np.zeros(0)
        self.tree = cKDTree(self.points) if len(self.points) else None
        self.segments = segments

    def query(self, points, max_dist=2.5):
        """(distance, bearing) of the nearest sampled segment point; distance is inf beyond max_dist."""
        points = np.asarray(points, float)
        if self.tree is None:
            return np.full(len(points), np.inf), np.zeros(len(points))
        d, j = self.tree.query(points, distance_upper_bound=max_dist)
        j = np.where(np.isfinite(d), j, 0)
        return d, self.bearings[j]
