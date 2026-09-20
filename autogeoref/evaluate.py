"""Leave-one-out evaluation against the team's own placements, and threshold calibration.

The team's files are affine and vertex-edited, so the truth pose is the rigid fit of the sheet to
the manual file and its tolerance is that fit's own rms, never tighter than 1.5 m / 1.5 degrees.
Parcels whose rigid rms exceeds REPORT_ONLY_RMS (47B, 48B, 42A, 43B, 169 on Kizhikaranai; 55 and 56
on Thailavaram) are reported but not counted: their own truth is too loose to judge anything by.

The error measured is the distance between the placed sheet and the sheet at the truth pose. It is
NOT the solver's `shift_m`, which only says how far the adjustment moved the starting guess.
"""
import csv
import datetime
import shutil

import numpy as np

from . import anchors, engine, fit, neighbours, paths, review, sheets

REPORT_ONLY_RMS = 3.0


def truth_pose(village, survey):
    pose = anchors.pose_from_points(village, survey)
    source = "points"
    if pose is None or pose[2] > REPORT_ONLY_RMS:
        hand = paths.manual_files(village, survey)
        if not hand:
            return None
        pose = anchors.pose_from_geometry(village, survey, hand[0])
        source = "geometry"
    if pose is None:
        return None
    theta, t, rms, scale, aniso, n = pose
    return {"theta": float(theta), "t": np.asarray(t, float), "rms": float(rms),
            "scale": float(scale), "anisotropy": float(aniso), "n": int(n),
            "heading_spread": _jackknife_heading(village, survey), "source": source}


def _jackknife_heading(village, survey):
    """How much the recovered heading moves when one control point is dropped."""
    p = paths.points_path(village, survey)
    if not p.exists():
        return 0.0
    P = anchors.read_points(p)
    if len(P) < 4:
        return 0.0
    out = []
    for k in range(len(P)):
        keep = np.delete(P, k, axis=0)
        theta, _t, _rms, _mx = fit.rigid_fit(keep[:, 2:4], keep[:, 0:2])
        out.append(theta)
    return float(max(out) - min(out))


def tolerances(truth):
    return max(1.5, float(truth["rms"])), max(1.5, float(truth.get("heading_spread", 0.0)))


def centroid_error(placed, truth_geom):
    return float(placed.centroid.distance(truth_geom.centroid))


def heading_error(got_deg, truth_deg):
    return float(abs(((got_deg - truth_deg) + 180) % 360 - 180))


def _sheet_at(village, survey, theta, t):
    return sheets.dissolve([(None, fit.apply_pose(g, theta, t))
                            for _p, g in sheets.load_sheet(village, survey)])


def leave_one_out(village, surveys=None, mode="full", keep_work=True):
    """Hide one manual parcel at a time on a copy of the village folder and score the placement."""
    stamp = datetime.datetime.now().strftime("%Y%m%d")
    work_root = paths.logs_dir() / ("loo_%s_%s" % (village, stamp))
    src = paths.vector_dir(village)
    real_project = paths.PROJECT
    targets = surveys or sorted(
        {p.name.split("_parcels_modified")[0] for p in src.glob("*_parcels_modified*.gpkg")},
        key=paths.survey_sort_key)
    results = []
    for survey in targets:
        truth = truth_pose(village, survey)
        if truth is None:
            continue
        truth_geom = _sheet_at(village, survey, truth["theta"], truth["t"])
        work = work_root / survey
        if work.exists():
            shutil.rmtree(work)
        (work / "FMB_Vector" / village).mkdir(parents=True)
        for f in src.iterdir():
            if f.is_dir():
                if f.name in ("neighbour_transcription", "gcp"):
                    shutil.copytree(f, work / "FMB_Vector" / village / f.name, dirs_exist_ok=True)
                continue
            if f.name.startswith(survey + "_parcels_modified"):
                continue                                   # the parcel under test loses its anchor
            if f.name == ("%s_parcels.geojson.points" % survey):
                continue                                   # and its control points
            if f.name in ("georef_status.csv", "anchors.csv"):
                continue                                   # a fresh run, not the last one's verdicts
            shutil.copy2(f, work / "FMB_Vector" / village / f.name)
        # the raster pin holds an absolute path, so copying the small json reuses the one GeoTIFF
        # instead of pulling 67 MB of Google tiles again for every survey under test
        georef_src = real_project / "FMB_Georef" / village
        if georef_src.exists():
            (work / "FMB_Georef" / village).mkdir(parents=True, exist_ok=True)
            for f in georef_src.glob("*"):
                if f.is_file() and f.suffix in (".json", ".xml"):
                    shutil.copy2(f, work / "FMB_Georef" / village / f.name)
        sigma0 = engine.SIGMA_AUTO
        try:
            paths.PROJECT = work
            neighbours._CACHE.pop(village, None)            # keyed by village, not by project root
            if mode == "gcp_only":
                engine.SIGMA_AUTO = 1e6                    # anchors contribute nothing
            engine.run(village, do_raster=False, do_topology=False, do_review=False)
            row = next((r for r in review.read_status(village) if r["survey"] == survey), None)
            placed = None
            out_gpkg = paths.output_path(village, survey)
            if out_gpkg.exists():
                import geopandas as gpd
                got = gpd.read_file(out_gpkg, layer="parcels")
                placed = sheets.dissolve([(None, g) for g in got.geometry])
        finally:
            paths.PROJECT = real_project
            engine.SIGMA_AUTO = sigma0
            neighbours._CACHE.pop(village, None)
        tol_m, tol_deg = tolerances(truth)
        base = {"survey": survey, "mode": mode, "truth_rms_m": round(truth["rms"], 2),
                "truth_source": truth["source"], "tolerance_m": round(tol_m, 2),
                "tolerance_deg": round(tol_deg, 2),
                "counted": truth["rms"] <= REPORT_ONLY_RMS}
        if row is None or placed is None:
            results.append(dict(base, colour="missing", centroid_error_m=None,
                                heading_error_deg=None, passed=False))
            continue
        err_m = centroid_error(placed, truth_geom)
        err_deg = heading_error(float(row.get("heading_deg") or 0.0), truth["theta"])
        results.append(dict(base, colour=row.get("colour"), method=row.get("method"),
                            share=row.get("share"), margin=row.get("margin"),
                            observable=row.get("observable"),
                            boundary_rms_m=row.get("boundary_rms_m"),
                            neighbours=row.get("neighbours"),
                            centroid_error_m=round(err_m, 2),
                            heading_error_deg=round(err_deg, 2),
                            passed=bool(err_m <= tol_m and err_deg <= tol_deg)))
        if not keep_work:
            shutil.rmtree(work, ignore_errors=True)
    return results


def ring_around(village, seeds, rings=1):
    """Surveys reachable from the seeds through the printed-neighbour lists in `rings` steps."""
    have = set(paths.surveys_with_sheets(village))
    frontier, reached = {str(s) for s in seeds}, set()
    for _ in range(rings):
        nxt = set()
        for s in frontier:
            for o in neighbours.printed(village, s):
                if o in have and o not in reached and o not in seeds:
                    nxt.add(o)
            for o in have:
                if s in neighbours.printed(village, o) and o not in reached and o not in seeds:
                    nxt.add(o)
        reached |= nxt
        frontier = nxt
    return sorted(reached, key=paths.survey_sort_key)


def seed_run(village, seeds, do_topology=True, label=None, rings=None, raster_margin_m=60.0):
    """Keep only `seeds` as hand-placed parcels on a copy of the village and place everything
    else from them. This is the rollout situation: a village with one or two seeds.

    Returns (rows, work_dir); each row has the tool's verdict plus the error against the team's
    own placement where one exists. Nothing under the real village folder is touched.
    """
    seeds = [str(s) for s in seeds]
    stamp = datetime.datetime.now().strftime("%Y%m%d")
    work = paths.logs_dir() / ("seed_%s_%s_%s" % (village, label or "-".join(seeds), stamp))
    src = paths.vector_dir(village)
    real_project = paths.PROJECT
    if work.exists():
        shutil.rmtree(work)
    (work / "FMB_Vector" / village).mkdir(parents=True)
    for f in src.iterdir():
        if f.is_dir():
            if f.name in ("neighbour_transcription",):
                shutil.copytree(f, work / "FMB_Vector" / village / f.name, dirs_exist_ok=True)
            continue
        name = f.name
        if "_parcels_modified" in name or name.endswith(".points"):
            survey = name.split("_parcels")[0]
            if survey not in seeds:
                continue                              # only the seeds keep their hand placement
        if name in ("georef_status.csv", "anchors.csv", "rail_conflicts.csv"):
            continue
        shutil.copy2(f, work / "FMB_Vector" / village / name)
    georef_src = real_project / "FMB_Georef" / village
    if georef_src.exists():
        (work / "FMB_Georef" / village).mkdir(parents=True, exist_ok=True)
        for f in georef_src.glob("*"):
            if f.is_file() and f.suffix in (".json", ".xml"):
                shutil.copy2(f, work / "FMB_Georef" / village / f.name)
    truths = {}
    for s in paths.surveys_with_sheets(village):
        if s in seeds:
            continue
        tr = truth_pose(village, s)
        if tr is not None:
            truths[s] = (tr, _sheet_at(village, s, tr["theta"], tr["t"]))
    only = None
    raster_bounds = None
    if rings is not None:
        only = ring_around(village, seeds, rings)
        # a satellite window just around the seed and its ring, not the whole village
        seed_bodies = []
        for s in seeds:
            tr = truth_pose(village, s)
            if tr is not None:
                seed_bodies.append(_sheet_at(village, s, tr["theta"], tr["t"]))
        if seed_bodies:
            from shapely.ops import unary_union
            b = unary_union(seed_bodies).bounds
            # neighbours the ring may reach lie within roughly one parcel width of the seed
            m = raster_margin_m + 0.5 * max(b[2] - b[0], b[3] - b[1])
            raster_bounds = (b[0] - m, b[1] - m, b[2] + m, b[3] + m)
    try:
        paths.PROJECT = work
        neighbours._CACHE.pop(village, None)
        if raster_bounds is not None:
            pin = work / "FMB_Georef" / village / "raster.json"
            if pin.exists():
                pin.unlink()                    # force a fresh, local export in the work copy
        summary = engine.run(village, do_raster=raster_bounds is not None, do_topology=do_topology,
                             do_review=False, only=only, raster_bounds=raster_bounds)
        summary["only"] = only
        rows = review.read_status(village)
        for r in rows:
            r["error_m"] = ""
            r["heading_error_deg"] = ""
            out = paths.output_path(village, r["survey"])
            if r["survey"] in truths and out.exists():
                import geopandas as gpd
                got = gpd.read_file(out, layer="parcels")
                placed = sheets.dissolve([(None, g) for g in got.geometry])
                tr, tg = truths[r["survey"]]
                r["error_m"] = round(centroid_error(placed, tg), 2)
                r["heading_error_deg"] = round(heading_error(float(r.get("heading_deg") or 0.0), tr["theta"]), 2)
    finally:
        paths.PROJECT = real_project
        neighbours._CACHE.pop(village, None)
    return rows, work, summary


def calibrate(rows):
    """Green must sit above every best-wrong score; amber above the median wrong score."""
    wrong = [r["share_best_wrong"] for r in rows if r.get("share_best_wrong") is not None]
    green = max(wrong) + 0.05 if wrong else 0.60
    amber = float(np.median(wrong)) if wrong else 0.35
    green = round(min(green, 0.95), 3)
    return {"green_share": green, "amber_share": round(min(amber, green - 0.05), 3),
            "table": sorted(rows, key=lambda r: -(r.get("share_truth") or 0))}


def write_report(village, rows, name="leave_one_out"):
    p = paths.logs_dir() / ("%s_%s_%s.csv" % (name, village,
                                              datetime.datetime.now().strftime("%Y%m%d")))
    p.parent.mkdir(parents=True, exist_ok=True)
    cols = sorted({k for r in rows for k in r})
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return p
