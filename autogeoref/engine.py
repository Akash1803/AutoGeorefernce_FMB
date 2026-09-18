"""The per-village pipeline of spec section 5.

Order: freeze the team's work, pick a start pose for everything else, let the satellite imagery
correct that pose where the evidence is unambiguous, tie the free parcels to their placed
neighbours in one block adjustment, write the files, then report. Anchors are never moved and
never rewritten; Puvi only ever suggests a starting guess and a reporting distance.
"""
import datetime
import math

import geopandas as gpd
import numpy as np
from pyproj import Proj
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from . import (align, anchors, edges as edgemod, files, fit, gcp, match, neighbours,
               paths, raster, review, sheets, topology, window)

ACRE = 0.000247105381
SIGMA_AUTO = 0.30
SIGMA_IMAGE = 1.0
SIGMA_START = 10.0
IMAGE_MAX_SHIFT = 1.5        # imagery may move an anchored pose by at most this
IMAGE_MAX_TURN = 2.0
GREEN_SHARE = 0.60           # recalibrated in Task 17 and written back here
AMBER_SHARE = 0.35


def colour_of(row):
    if row.get("method") == "puvi-only" or row.get("ambiguous"):
        return "red"
    rms = row.get("boundary_rms_m")
    allow = SIGMA_AUTO + (row.get("anchor_rms_m") or 0.0)
    anchored = rms is not None and rms <= allow
    strong_image = (row.get("share") or 0.0) >= GREEN_SHARE and row.get("observable")
    if anchored and (strong_image or (row.get("n_neighbours") or 0) >= 2):
        return "green"
    if anchored or strong_image or (row.get("share") or 0.0) >= AMBER_SHARE:
        return "amber"
    return "red"


def write_parcels(village, survey, theta, t, row, out_dir=None):
    """<survey>_parcels_modified.gpkg with parcels + edges layers and the section 11 attributes."""
    sheet = sheets.load_sheet(village, survey)
    proj = Proj("EPSG:32644")
    run = row.get("run") or datetime.datetime.now().isoformat(timespec="seconds")
    parcels, edge_rows = [], []
    for props, geom in sheet:
        placed = fit.apply_pose(geom, theta, t)
        fit.assert_rigid(geom, placed)
        area = placed.area
        rec = dict(props)
        rec.update({"area_sqm": round(area, 3), "area_are": round(area / 100, 4),
                    "area_hect": round(area / 1e4, 6), "area_acre": round(area * ACRE, 5),
                    "area_cent": round(area * ACRE * 100, 3), "perimeter_m": round(placed.length, 3),
                    "fmb_area_sqm": round(geom.area, 3), "fmb_area_acre": round(geom.area * ACRE, 5),
                    "fmb_perimeter_m": round(geom.length, 3),
                    "georef_method": row.get("method"), "georef_colour": row.get("colour"),
                    "georef_confidence": row.get("confidence"),
                    "georef_residual_m": row.get("boundary_rms_m"),
                    "anchored_to": row.get("neighbours"), "heading_deg": round(theta, 3),
                    "puvi_reference_m": row.get("puvi_reference_m"), "georef_notes": row.get("notes"),
                    "georef_run": run})
        parcels.append({**rec, "geometry": placed})
        ring = list(placed.exterior.coords)
        if Polygon(ring).exterior.is_ccw:
            ring = ring[::-1]
        ring = ring[:-1]
        start = max(range(len(ring)), key=lambda i: ring[i][1])
        ring = ring[start:] + ring[:start]
        for i in range(len(ring)):
            (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % len(ring)]
            grid = math.degrees(math.atan2(x1 - x0, y1 - y0)) % 360
            lon, lat = proj((x0 + x1) / 2, (y0 + y1) / 2, inverse=True)
            conv = proj.get_factors(lon, lat).meridian_convergence
            edge_rows.append({"survey_no": survey, "poly_id": props.get("poly_id"),
                              "plot_no": props.get("plot_no"), "edge_no": i + 1,
                              "length_m": round(math.dist((x0, y0), (x1, y1)), 3),
                              "bearing_grid_deg": round(grid, 2),
                              "bearing_true_deg": round((grid + conv) % 360, 2),
                              "grid_convergence_deg": round(conv, 3),
                              "geometry": LineString([(x0, y0), (x1, y1)])})
    target = ((out_dir / ("%s_parcels_modified.gpkg" % survey)) if out_dir
              else paths.output_path(village, survey))

    def writer(path):
        gpd.GeoDataFrame(parcels, geometry="geometry", crs="EPSG:32644").to_file(
            path, layer="parcels", driver="GPKG")
        gpd.GeoDataFrame(edge_rows, geometry="geometry", crs="EPSG:32644").to_file(
            path, layer="edges", driver="GPKG")

    if target.exists():
        if not files.safe_write(target, writer):
            raise PermissionError("cannot replace %s (open in QGIS?)" % target)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        writer(target)
    return target


def run(village, do_raster=False, do_topology=False, do_review=True, project=None):
    """Steps 1-13 of spec section 5. Returns a summary dict; details go to georef_status.csv."""
    run_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    tool_fps = review.tool_written_fingerprints(village)
    anchor_map = anchors.load_anchors(village, tool_fps)
    conflicts = anchors.anchor_conflicts(village, anchor_map)
    superseded = {s: paths.manual_files(village, s)[1:] for s in anchor_map}
    anchors.write_anchors_csv(village, anchor_map, superseded, conflicts)

    todo = [s for s in paths.surveys_with_sheets(village) if s not in anchor_map]
    placed = {s: (a.theta, a.t) for s, a in anchor_map.items()}
    rows = []

    if todo and (do_raster or raster.pinned(village) is None):
        geoms = [anchors.placed_geometry(village, a) for a in anchor_map.values()]
        geoms = [g for g in geoms if g is not None and not g.is_empty]
        if geoms:
            b = unary_union(geoms).bounds
            raster.export(village, (b[0] - 100, b[1] - 100, b[2] + 100, b[3] + 100), run_id=run_id)

    index = None
    if todo and raster.pinned(village) is not None:
        index = edgemod.SegmentIndex(edgemod.detect(raster.pinned(village)["path"], min_len_m=3.0))

    for survey in todo:
        row = {"survey": survey, "run": stamp, "method": "puvi-only", "ambiguous": True,
               "notes": "", "n_neighbours": 0, "file": paths.output_path(village, survey).name}
        hyps, _win = window.start_poses(village, survey)
        pose = hyps[0] if hyps else None
        if index is not None and hyps:
            pts, brg = sheets.sample_outline(sheets.outline(sheets.load_sheet(village, survey)), 1.0)
            res = align.search(pts, brg, index, hyps)
            row.update({"share": round(res.share, 3), "margin": round(res.margin, 3),
                        "observable": res.observable, "ambiguous": res.ambiguous, "notes": res.note})
            if not res.ambiguous and res.observable:
                pose = (res.theta, res.t)
                row["method"] = "image"
        if pose is None:
            rows.append(dict(row, colour="red", confidence=0, status="waiting",
                             notes=(row["notes"] + " | no Puvi polygon and no placed neighbour").strip(" |")))
            continue
        row["heading_deg"] = round(pose[0], 3)
        placed[survey] = pose
        rows.append(row)

    # block adjustment: anchors fixed, everything else free
    free = {s: placed[s] for s in todo if s in placed}
    pair_obs, line_obs, gcp_obs, priors = [], [], [], []
    for s in free:
        priors.append(fit.PosePrior(s, free[s][0], tuple(free[s][1]), SIGMA_START, 10.0))
        v_s = sheets.outline(sheets.load_sheet(village, s))
        v_s_ground = [tuple(p) for p in fit.transform_points(v_s, free[s][0], free[s][1])]
        for other, (th_o, t_o) in placed.items():
            if other == s or neighbours.are_neighbours(village, s, other) is False:
                continue
            v_o = [tuple(p) for p in fit.transform_points(
                sheets.outline(sheets.load_sheet(village, other)), th_o, t_o)]
            for chain in match.common_chains(v_s_ground, v_o)[:1]:
                sigma = SIGMA_AUTO + (anchor_map[other].rms if other in anchor_map else 0.0)
                pair_obs += match.chain_observations(chain, s, other, sigma)
    adjusted = fit.block_adjust(free, {s: (a.theta, a.t) for s, a in anchor_map.items()},
                                pair_obs, line_obs, gcp_obs, priors)
    poses = {s: {"theta": v["theta"], "t": v["t"]} for s, v in adjusted.items()}
    poses.update({s: {"theta": a.theta, "t": a.t} for s, a in anchor_map.items()})
    residuals = fit.residuals_by_pair(poses, pair_obs)

    written = []
    for row in rows:
        s = row["survey"]
        if s not in adjusted:
            continue
        theta, t = adjusted[s]["theta"], adjusted[s]["t"]
        mine = [v for k, v in residuals.items() if s in k]
        row["boundary_rms_m"] = round(min(v["rms"] for v in mine), 3) if mine else None
        partners = sorted({k[0] if k[1] == s else k[1] for k in residuals if s in k},
                          key=paths.survey_sort_key)
        row["n_neighbours"] = len(partners)
        row["neighbours"] = ",".join(partners)
        row["anchor_rms_m"] = max([anchor_map[n].rms for n in partners if n in anchor_map] or [0.0])
        row["shift_m"] = round(adjusted[s]["shift_m"], 2)
        row["sigma_pos_m"] = round(adjusted[s]["sigma_pos_m"], 3)
        row["sigma_head_deg"] = round(adjusted[s]["sigma_head_deg"], 3)
        if row["boundary_rms_m"] is not None and row["method"] != "image":
            row["method"] = "anchors"
        body = unary_union([fit.apply_pose(g, theta, t) for _p, g in sheets.load_sheet(village, s)])
        row["puvi_reference_m"] = round(window.puvi_distance(village, s, body) or 0.0, 1)
        trusted, ratio = window.puvi_trusted(village, s, sheets.dissolve(sheets.load_sheet(village, s)).area)
        row["puvi_trusted"] = trusted
        row["puvi_ratio"] = round(ratio, 2) if ratio == ratio else ""
        row["colour"] = colour_of(row)
        row["confidence"] = {"green": 90, "amber": 60, "red": 30}[row["colour"]]
        row["status"] = "placed"
        target = write_parcels(village, s, theta, t, row)
        row["fp_placed"] = anchors.fingerprint(target)
        if index is not None:
            gcp.write(village, s, gcp.corner_gcps(village, s, theta, t, index), _crs_wkt())
        written.append(s)

    if do_topology and written:
        geoms = {s: unary_union([fit.apply_pose(g, adjusted[s]["theta"], adjusted[s]["t"])
                                 for _p, g in sheets.load_sheet(village, s)]) for s in written}
        geoms.update({s: anchors.placed_geometry(village, a) for s, a in anchor_map.items()})
        rail = {s for s in geoms if s in _rail_parcels(village)}
        _fixed, report = topology.fix(geoms, movable=set(written), rail=rail)
        topology.report_rail_conflicts(village, report)

    for row in rows:
        if row.get("fp_placed"):
            row["fp_final"] = row["fp_placed"]
    review.write_status(village, rows)
    if do_review:
        review.build_group(village, rows, project)
    return {"village": village, "anchors": len(anchor_map), "placed": len(written),
            "waiting": sum(1 for r in rows if r.get("status") == "waiting"),
            "conflicts": len(conflicts), "run": run_id, "rows": rows}


def _crs_wkt():
    from pyproj import CRS
    return CRS.from_epsg(32644).to_wkt()


def _rail_parcels(village):
    """Surveys whose buffer polygon carries the traced rail centreline (reference knowledge)."""
    try:
        buf = gpd.read_file(paths.BUFFER_GPKG, layer=paths.BUFFER_LAYER).to_crs(32644)
        rail = unary_union(list(gpd.read_file(paths.RAIL_GPKG, layer=paths.RAIL_LAYER).to_crs(32644).geometry))
    except Exception:
        return set()
    sub = buf[buf["village_code"].astype(str).str.strip() == village]
    return {str(r.survey_no).strip() for r in sub.itertuples()
            if r.geometry.intersection(rail).length >= 50.0}
