"""Fit a sheet outline to the straight edges visible in the satellite raster.

Two guards make this safe on a railway corridor, where the 2026-09-18 trial showed a free image
search sliding up to 28 m or flipping 180 degrees:

* margin - the best pose must beat the best pose more than EXCL_M / EXCL_DEG away by MARGIN_MIN,
  and no local maximum within PEAK_RATIO of the best may sit more than 1 m from it;
* observability - the matched edges must span two directions (at least CROSS_ANGLE apart),
  otherwise the along-track position is not observable from imagery and a neighbour must supply it.

The score is the share of 1 m outline samples lying within MATCH_DIST of a same-direction segment,
which is the metric the spike measured, so calibrated thresholds apply to it directly.
"""
import math
from dataclasses import dataclass

import numpy as np

from . import fit

MATCH_DIST = 1.5          # metres: the spike's own tolerance
MATCH_DEG = 12.0          # degrees
EXCL_M = 3.0              # margin exclusion zone
EXCL_DEG = 3.0
MARGIN_MIN = 0.15         # the best pose must beat the best outside pose by this much share
PEAK_RATIO = 0.75         # local maxima above this fraction of the best count as competitors
PRIMARY_MIN = 40.0        # metres of matched outline in the main direction
CROSS_MIN = 20.0          # metres matched in a direction at least CROSS_ANGLE from it
CROSS_ANGLE = 30.0        # degrees
GROUP_TOL = 10.0          # degrees: how close two sample bearings must be to count as one direction


@dataclass
class AlignResult:
    theta: float
    t: np.ndarray
    share: float
    margin: float
    alternative: tuple
    ambiguous: bool
    matched_by_direction: list
    observable: bool
    note: str


def _matched_mask(samples, bearings, index, theta_deg, t):
    pts = fit.transform_points(samples, theta_deg, t)
    d, b = index.query(pts, max_dist=MATCH_DIST + 1.0)
    turned = (bearings + math.radians(theta_deg)) % math.pi
    diff = np.abs((b - turned + math.pi / 2) % math.pi - math.pi / 2)
    return np.isfinite(d) & (d <= MATCH_DIST) & (diff <= math.radians(MATCH_DEG))


def share(samples, bearings, index, theta_deg, t):
    """Fraction of outline samples lying on a same-direction image edge."""
    if len(samples) == 0:
        return 0.0
    return float(_matched_mask(samples, bearings, index, theta_deg, t).mean())


def observability(samples, bearings, index, theta_deg, t):
    """How much matched outline length lies in the main direction and across it.

    Samples are one metre apart, so counting them counts metres.
    """
    mask = _matched_mask(samples, bearings, index, theta_deg, t)
    turned = ((bearings + math.radians(theta_deg)) % math.pi)[mask]
    if turned.size == 0:
        return {"primary_m": 0.0, "cross_m": 0.0, "cross_angle_deg": 0.0, "ok": False}
    groups = []
    for b in turned:
        for i, (gb, n) in enumerate(groups):
            if min(abs(b - gb), math.pi - abs(b - gb)) <= math.radians(GROUP_TOL):
                groups[i] = (gb, n + 1)
                break
        else:
            groups.append((b, 1))
    groups.sort(key=lambda g: -g[1])
    primary_b, primary_n = groups[0]
    cross = [(gb, n) for gb, n in groups[1:]
             if min(abs(gb - primary_b), math.pi - abs(gb - primary_b)) >= math.radians(CROSS_ANGLE)]
    cross_b, cross_n = (cross[0] if cross else (primary_b, 0))
    ang = math.degrees(min(abs(cross_b - primary_b), math.pi - abs(cross_b - primary_b)))
    return {"primary_m": float(primary_n), "cross_m": float(cross_n), "cross_angle_deg": round(ang, 1),
            "ok": bool(primary_n >= PRIMARY_MIN and cross_n >= CROSS_MIN)}


def search(samples, bearings, index, hypotheses, shift_m=10.0, rot_deg=10.0,
           coarse_m=1.0, coarse_deg=1.0, fine_m=0.25, fine_deg=0.25):
    """Best pose around the given hypotheses, with the margin and observability verdicts."""
    grid = []
    for theta0, t0 in hypotheses:
        t0 = np.asarray(t0, float)
        for dth in np.arange(-rot_deg, rot_deg + 1e-9, coarse_deg):
            for dx in np.arange(-shift_m, shift_m + 1e-9, coarse_m):
                for dy in np.arange(-shift_m, shift_m + 1e-9, coarse_m):
                    th = theta0 + dth
                    t = t0 + np.array([dx, dy])
                    grid.append((share(samples, bearings, index, th, t), th, t))
    if not grid:
        return AlignResult(0.0, np.zeros(2), 0.0, 0.0, None, True, [], False, "no hypothesis")
    grid.sort(key=lambda g: -g[0])
    best_s, best_th, best_t = grid[0]
    for dth in np.arange(-coarse_deg, coarse_deg + 1e-9, fine_deg):
        for dx in np.arange(-coarse_m, coarse_m + 1e-9, fine_m):
            for dy in np.arange(-coarse_m, coarse_m + 1e-9, fine_m):
                th = best_th + dth
                t = best_t + np.array([dx, dy])
                s = share(samples, bearings, index, th, t)
                if s > best_s:
                    best_s, best_th, best_t = s, th, t
    competitor = None
    for s, th, t in grid:
        if (float(np.linalg.norm(t - best_t)) <= EXCL_M
                and abs(((th - best_th) + 180) % 360 - 180) <= EXCL_DEG):
            continue
        competitor = (s, th, t)
        break
    margin = best_s - (competitor[0] if competitor else 0.0)
    # a competing peak, not the shoulder of the single peak we are on: the score surface is about
    # MATCH_DIST wide, so only high scores beyond the exclusion zone count as a rival
    near_peaks = [(s, th, t) for s, th, t in grid
                  if s >= PEAK_RATIO * best_s and float(np.linalg.norm(t - best_t)) > EXCL_M]
    ambiguous = bool(margin < MARGIN_MIN or near_peaks)
    obs = observability(samples, bearings, index, best_th, best_t)
    note = ""
    if ambiguous and competitor is not None:
        dx, dy = competitor[2] - best_t
        note = "alternative pose %.1f m away at bearing %.0f deg" % (
            math.hypot(dx, dy), math.degrees(math.atan2(dx, dy)) % 360)
    elif not obs["ok"]:
        note = ("matched edges in one direction only (%.0f m); along-track from neighbours"
                % obs["primary_m"])
    return AlignResult(float(best_th), np.asarray(best_t, float), float(best_s), float(margin),
                       (competitor[1], competitor[2]) if competitor else None, ambiguous,
                       [obs["primary_m"], obs["cross_m"], obs["cross_angle_deg"]], obs["ok"], note)
