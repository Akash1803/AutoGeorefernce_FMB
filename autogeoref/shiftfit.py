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

Only a translation is applied. Shape, size, rotation and area are left exactly as Puvi drew them,
which is what the task asked for. Nothing existing is deleted or edited: outputs are new files.
"""
import csv
import datetime
import logging

import numpy as np

from . import paths

log = logging.getLogger(__name__)

IDW_MIN_CONTROL = 4     # below this a local field is guesswork: one mean shift is used instead
IDW_K = 4               # control points that vote for a target
IDW_POWER = 2.0         # inverse-distance weight exponent
FAR_CONTROL_M = 1500.0  # a control point further than this from a target only votes if nothing closer


def normalise(survey):
    """Survey keys across the three sources: Puvi writes '10 B', the portal and sheets '10B'."""
    return str(survey).strip().replace(" ", "").upper()


def mean_shift(disp):
    """The one translation that best explains every control displacement."""
    return np.asarray(disp, float).mean(axis=0)


def idw(points, disp, target, k=IDW_K, power=IDW_POWER):
    """Inverse-distance weighted displacement at `target`, from control at `points`."""
    points = np.asarray(points, float)
    disp = np.asarray(disp, float)
    d = np.linalg.norm(points - np.asarray(target, float), axis=1)
    if len(points) == 0:
        return np.zeros(2)
    if (d < 1e-6).any():
        return disp[int(np.argmin(d))]
    near = np.argsort(d)[:min(k, len(points))]
    w = 1.0 / d[near] ** power
    return (w @ disp[near]) / w.sum()


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


REPORT_COLUMNS = ["village_code", "village_name", "control", "targets", "method",
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
