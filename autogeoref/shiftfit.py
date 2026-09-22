"""Shift the Puvi cadastral vectors onto the ground using the team's own placements as control.

Why this and not the sheet engine: the engine places an FMB sheet by chaining it off a parcel that
is already on the ground, which needs a seed and a transcription per village. This module answers a
smaller question, the one asked on 2026-09-22: move the Puvi survey polygons so they sit closer to
where the parcels really are, quickly, for the whole rail buffer, without touching their shape.

Measured on 2026-09-22 against the team's hand placements (median error of a corrected parcel):

| village      | Puvi as it is | 1 control | 2 | 3 | 5 | 10 | all (leave-one-out) |
|--------------|---------------|-----------|---|---|---|----|---------------------|
| Thailavaram  | 6.4 m         | 6.4       | 5.2 | 5.0 | 4.6 | 4.1 | 3.4 |
| Kizhikaranai | 8.3 m         | 8.2       | 6.9 | 7.5 | 6.5 | 6.8 | 5.5 |
| Thirukatchur | 24.1 m        | 10.3      | 11.4 | 9.6 | 6.2 | 5.0 | 4.0 |

Three automatic alternatives were tried first and all failed, so none of them is in here: matching a
Puvi polygon to satellite edges (median error 8.3 m -> 24.9 m, every result ambiguous), matching a
whole village block (share 0.02, drifts to the edge of the search box), and matching an accurate
sheet outline positioned by Puvi (8.3 m -> 17.6 m, margin under 0.12 everywhere). The satellite can
confirm a placement; at this image quality it cannot decide one within 20 m.

Parcels are moved through one smooth field rather than each by its own vector. Moving each parcel
on its own splits every boundary two of them share, which is the double line Akash saw on
2026-09-22: 55 overlapping pairs in Thirukatchur where Puvi had none, and the length of shared
drawn lines halved. Evaluating the field at each distinct vertex instead keeps the fabric closed:
zero overlaps, zero new gaps, shared lines back to 9146 m of Puvi's 9095 m. The cost is that a
parcel is no longer moved perfectly rigidly; with the tuned field the area change is 2 % at the
90th percentile and 8.3 % at worst, against 9 to 18 % and 57 % before tuning.

Nothing existing is deleted or edited: outputs are new files.
"""
import csv
import datetime
import logging

import numpy as np

from . import paths

log = logging.getLogger(__name__)

IDW_MIN_CONTROL = 4     # below this a local field is guesswork: one mean shift is used instead
IDW_K = 12              # control points that vote for a point
IDW_POWER = 1.0         # inverse-distance weight exponent
IDW_SMOOTH_M = 60.0     # a control point never outvotes the rest from closer than this

# Why these and not a sharper field: a sharp field pulls the near side of a parcel much harder than
# the far side, which distorts it. Measured on 2026-09-22 over the three villages that have control
# (leave-one-out error, then the 90th percentile of the area change the warp causes):
#   k=4  power=2 no smoothing : 4.0 / 3.6 / 5.5 m, area 9.0 / 6.8 / 18.2 %, worst 57 %
#   k=12 power=1 smooth 60 m  : 3.8 / 3.6 / 6.0 m, area 2.2 / 1.7 /  2.2 %, worst  8 %
#   all  power=1 smooth 120 m : 5.1 / 4.2 / 5.2 m, area 1.4 / 0.6 /  1.2 %, worst  4 %
# The middle row is as accurate as the sharp field and holds the shapes, so it is the default.


def normalise(survey):
    """Survey keys across the three sources: Puvi writes '10 B', the portal and sheets '10B'."""
    return str(survey).strip().replace(" ", "").upper()


def mean_shift(disp):
    """The one translation that best explains every control displacement."""
    return np.asarray(disp, float).mean(axis=0)


def idw(points, disp, target, k=IDW_K, power=IDW_POWER, smooth=IDW_SMOOTH_M):
    """Inverse-distance weighted displacement at `target`, from control at `points`.

    `smooth` is a floor on the distance, so a control point a metre away does not drown out the
    twelve around it: without it the field has a spike at every control parcel and a parcel sitting
    on one is stretched (up to 57 % of its area on the first run).
    """
    points = np.asarray(points, float)
    disp = np.asarray(disp, float)
    if len(points) == 0:
        return np.zeros(2)
    d = np.linalg.norm(points - np.asarray(target, float), axis=1)
    near = np.argsort(d)[:min(k, len(points))]
    w = 1.0 / (d[near] ** power + smooth ** power)
    return (w @ disp[near]) / w.sum()


def warp_points(points, disp, coords, min_control=IDW_MIN_CONTROL):
    """Displacement for arbitrary coordinates, evaluated point by point.

    Moving whole parcels by one vector each splits every boundary two parcels share: the 2026-09-22
    run left 55 overlapping pairs in Thirukatchur where Puvi had none, and halved the length of
    shared drawn lines. Evaluating the same field at each vertex instead means two parcels that
    share a vertex share its movement, so the fabric stays closed.
    """
    return fit(points, disp, coords, min_control=min_control)


def fit(points, disp, targets, min_control=IDW_MIN_CONTROL):
    """Displacement for each target. Returns (array of shifts, method name).

    With few control points a local field would invent structure it cannot know, so one mean shift
    is applied to the whole village; with enough of them each target follows its nearest control.
    """
    points = np.asarray(points, float)
    disp = np.asarray(disp, float)
    targets = np.asarray(targets, float)
    if len(points) == 0:
        return np.zeros((len(targets), 2)), "none"
    if len(points) < min_control:
        return np.repeat(mean_shift(disp)[None, :], len(targets), axis=0), "mean shift"
    return np.array([idw(points, disp, t) for t in targets]), "local field"


def leave_one_out(points, disp, min_control=IDW_MIN_CONTROL):
    """Error the fit makes on control it has not seen. Returns the per-point error in metres."""
    points = np.asarray(points, float)
    disp = np.asarray(disp, float)
    if len(points) < 2:
        return np.array([])
    out = []
    for i in range(len(points)):
        rest = np.delete(np.arange(len(points)), i)
        pred, _ = fit(points[rest], disp[rest], points[i][None, :], min_control=min_control)
        out.append(float(np.linalg.norm(pred[0] - disp[i])))
    return np.array(out)


def accuracy(points, disp, min_control=IDW_MIN_CONTROL):
    """Summary of the leave-one-out error, and of the error the parcels started with."""
    loo = leave_one_out(points, disp, min_control=min_control)
    before = np.linalg.norm(np.asarray(disp, float), axis=1) if len(disp) else np.array([])
    return {
        "control": int(len(points)),
        "before_median_m": round(float(np.median(before)), 2) if len(before) else None,
        "before_p90_m": round(float(np.percentile(before, 90)), 2) if len(before) else None,
        "after_median_m": round(float(np.median(loo)), 2) if len(loo) else None,
        "after_p90_m": round(float(np.percentile(loo, 90)), 2) if len(loo) else None,
        "after_max_m": round(float(loo.max()), 2) if len(loo) else None,
    }


REPORT_COLUMNS = ["village_code", "village_name", "control", "targets", "method", "seam_note",
                  "before_median_m", "before_p90_m", "after_median_m", "after_p90_m", "after_max_m",
                  "mean_shift_x_m", "mean_shift_y_m", "run_at"]


def write_report(rows, path):
    """One row per village: how much control it had, and what the fit is worth."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in REPORT_COLUMNS})
    log.info("shift report -> %s", path)
    return path


def run_stamp():
    return datetime.datetime.now().isoformat(timespec="seconds")


def output_dir(stamp=None):
    """Where a corrected set is written: beside the project, never over the source."""
    stamp = stamp or datetime.date.today().strftime("%Y%m%d")
    return paths.PROJECT / "Puvi_Shifted" / stamp
