"""Puvi-derived search window and starting poses.

Puvi is a location REFERENCE ONLY (user rule): it may say roughly where a survey is and roughly
how it is turned, and its distance is reported. It never accepts, rejects or moves a placement.
Measured offsets on this corridor: 1.6 m to 36.8 m.

It is also not always the same parcel. For 22 corridor surveys the Puvi polygon and the FMB sheet
describe different land (Thirukatchur 569B: a 0.98 acre sheet against the 259 acre village tank),
so an area gate drops Puvi for that survey rather than seeding the search 700 m away.
"""
import math
import re

import geopandas as gpd
import numpy as np
from shapely import affinity
from shapely.ops import unary_union

from . import paths, sheets

BUFFER_M = 40.0        # railway strips were 37 m off; the window must still contain the truth
COARSE_STEP = 2.0      # degrees
FINE_STEP = 0.5        # degrees
SHIFT_STEP = 2.0       # metres
SHIFT_RANGE = 6.0      # metres
TWIN_RATIO = 1.15      # keep the 180 degree twin when it scores within 15 % of the best
AREA_GATE = 3.0        # Puvi and the portal describe different parcels beyond this factor
_CACHE = {}


def _puvi(village):
    if village not in _CACHE:
        shp = paths.puvi_shapefile(village)
        if shp is None:
            _CACHE[village] = None
        else:
            gdf = gpd.read_file(shp).to_crs(32644)
            gdf["key"] = gdf["survey_no"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
            _CACHE[village] = {k: unary_union(list(g.geometry)).buffer(0) for k, g in gdf.groupby("key")}
    return _CACHE[village]


def puvi_polygon(village, survey):
    """(geometry, key_kind) where key_kind is 'exact', 'parent' or ''."""
    table = _puvi(village)
    if not table:
        return None, ""
    key = str(survey).upper()
    if key in table:
        return table[key], "exact"
    m = re.match(r"^\d+", key)
    if m and m.group() in table:
        return table[m.group()], "parent"
    return None, ""


def puvi_trusted(village, survey, sheet_area):
    """(trusted, puvi_area / sheet_area). 22 corridor surveys fail this; 569B by a factor of 263."""
    target, _kind = puvi_polygon(village, survey)
    if target is None or sheet_area <= 0:
        return False, float("nan")
    ratio = target.area / sheet_area
    return bool((1.0 / AREA_GATE) <= ratio <= AREA_GATE), ratio


def start_poses(village, survey, buffer_m=BUFFER_M):
    """Hypotheses [(theta_deg, t)] from aligning the sheet outline to the Puvi shape, plus the window."""
    target, _kind = puvi_polygon(village, survey)
    if target is None:
        return [], None
    body = sheets.dissolve(sheets.load_sheet(village, survey))
    trusted, _ratio = puvi_trusted(village, survey, body.area)
    if not trusted:
        return [], None                 # different parcel: fall back to the neighbour window
    c = np.array(body.centroid.coords[0])
    tc = np.array(target.centroid.coords[0])

    def score(theta, dx=0.0, dy=0.0):
        g = affinity.rotate(body, theta, origin=(c[0], c[1]))
        g = affinity.translate(g, tc[0] - c[0] + dx, tc[1] - c[1] + dy)
        return g.symmetric_difference(target).area

    coarse = sorted(((score(th), th) for th in np.arange(0.0, 360.0, COARSE_STEP)), key=lambda x: x[0])
    best_s, best_th = coarse[0]
    twin_th = (best_th + 180.0) % 360.0
    twin_s = min(s for s, th in coarse if abs(((th - twin_th) + 180) % 360 - 180) <= COARSE_STEP)
    starts = [best_th] + ([twin_th] if twin_s <= best_s * TWIN_RATIO else [])

    poses = []
    for th0 in starts:
        best = (float("inf"), th0, 0.0, 0.0)
        for th in np.arange(th0 - 3.0, th0 + 3.0001, FINE_STEP):
            for dx in np.arange(-SHIFT_RANGE, SHIFT_RANGE + 0.001, SHIFT_STEP):
                for dy in np.arange(-SHIFT_RANGE, SHIFT_RANGE + 0.001, SHIFT_STEP):
                    s = score(th, dx, dy)
                    if s < best[0]:
                        best = (s, th, dx, dy)
        _s, th, dx, dy = best
        rad = math.radians(th)
        R = np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])
        poses.append((float(th), tc + np.array([dx, dy]) - R @ c))
    return poses, target.buffer(buffer_m)


def puvi_distance(village, survey, placed_geom):
    """Reporting only: how far the placement sits from the Puvi polygon of the same number."""
    target, _kind = puvi_polygon(village, survey)
    if target is None or placed_geom is None or placed_geom.is_empty:
        return None
    return float(placed_geom.centroid.distance(target.centroid))
