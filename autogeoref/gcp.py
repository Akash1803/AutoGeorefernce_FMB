"""Ground control points in the team's own format.

The team georeferences in the QGIS Georeferencer, which reads and writes
`<survey>_parcels.geojson.points`. The tool writes its automatic GCPs into that same file so the
team opens a sheet with the points already placed, drags the ones that are wrong, and saves.
Two dragged points fully determine a rigid pose, so there is no minimum-count rule for a refit.

A point is written for every real corner of the sheet, matched or not: a parcel the imagery could
not read is exactly the one the team needs points for.
"""
import csv
import datetime
from dataclasses import dataclass

import numpy as np

from . import anchors, fit, paths, sheets

CORNER_TURN_MIN = 20.0        # degrees: what counts as a real corner
MATCH_DIST = 1.5              # metres
TEAM_SIGMA = 0.20             # metres, spec section 5 step 11
CSV_COLUMNS = ["survey", "corner_id", "sheet_x", "sheet_y", "map_x", "map_y",
               "matched", "residual_m", "source", "written"]


@dataclass
class Gcp:
    survey: str
    corner_id: int
    sheet_xy: tuple
    map_xy: tuple
    matched: bool
    residual_m: float
    source: str               # "auto" | "team"


def corner_gcps(village, survey, theta, t, index, turn_min=CORNER_TURN_MIN):
    """One GCP per real corner of the sheet outline at the given pose, matched flag from imagery."""
    v = sheets.outline(sheets.load_sheet(village, survey))
    turns = sheets.turn_angles(v)
    out = []
    corner_id = 0
    for i, p in enumerate(v):
        if turns[i] < turn_min:
            continue
        corner_id += 1
        ground = fit.transform_points([p], theta, t)[0]
        d, _b = index.query(np.array([ground]), max_dist=MATCH_DIST + 1.0)
        matched = bool(np.isfinite(d[0]) and d[0] <= MATCH_DIST)
        out.append(Gcp(survey, corner_id, (float(p[0]), float(p[1])),
                       (float(ground[0]), float(ground[1])), matched,
                       float(d[0]) if np.isfinite(d[0]) else float("nan"), "auto"))
    return out


def _csv_rows(village, survey):
    p = paths.gcp_csv(village, survey)
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write(village, survey, gcps, crs_wkt):
    """Write the durable CSV always; write the .points file only when the team owns no version."""
    p = paths.gcp_csv(village, survey)
    p.parent.mkdir(parents=True, exist_ok=True)
    existing = {int(r["corner_id"]): r for r in _csv_rows(village, survey) if r.get("source") == "team"}
    rows = []
    for g in gcps:
        if g.corner_id in existing:          # never overwrite a team row
            rows.append(existing[g.corner_id])
            continue
        rows.append({"survey": g.survey, "corner_id": g.corner_id,
                     "sheet_x": "%.4f" % g.sheet_xy[0], "sheet_y": "%.4f" % g.sheet_xy[1],
                     "map_x": "%.4f" % g.map_xy[0], "map_y": "%.4f" % g.map_xy[1],
                     "matched": int(bool(g.matched)),
                     "residual_m": "" if g.residual_m != g.residual_m else "%.3f" % g.residual_m,
                     "source": g.source,
                     "written": datetime.datetime.now().isoformat(timespec="seconds")})
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    if team_points(village, survey) is None:
        anchors.write_points(paths.points_path(village, survey),
                             [(float(r["map_x"]), float(r["map_y"]),
                               float(r["sheet_x"]), float(r["sheet_y"])) for r in rows], crs_wkt)
    return p


def team_points(village, survey):
    """The team's GCPs when their .points file is newer than the tool's record, else None."""
    pts = paths.points_path(village, survey)
    if not pts.exists():
        return None
    record = paths.gcp_csv(village, survey)
    if record.exists() and pts.stat().st_mtime <= record.stat().st_mtime + 1:
        return None
    P = anchors.read_points(pts)
    if len(P) < 2:
        return None
    return [Gcp(survey, i + 1, (float(r[2]), float(r[3])), (float(r[0]), float(r[1])),
                True, float("nan"), "team") for i, r in enumerate(P)]


def refit(village, survey):
    """(theta, t, rms, n) from the team's GCPs, scale fixed at 1; None when there are fewer than two."""
    team = team_points(village, survey)
    if not team:
        return None
    P = np.array([g.sheet_xy for g in team], float)
    Q = np.array([g.map_xy for g in team], float)
    theta, t, rms, _mx = fit.rigid_fit(P, Q)
    return float(theta), t, float(rms), len(team)
