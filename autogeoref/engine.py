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
CHAIN_GATE_M = 40.0          # a chain whose paired points are further apart than the search window
                             # is a congruent boundary somewhere else, not this parcel's boundary
MAX_RESIDUAL_M = 5.0         # above this the placement contradicts its own neighbours
CHAIN_MIN_M = 12.0           # a shared run shorter than this pins nothing
CHAIN_FIT_RMS = 1.0          # metres: how well a chain must fit before its pose is a candidate
MAX_OVERLAP_SHARE = 0.02     # a pose may not sit on top of a parcel already on the ground
POSE_TOL_DEG = 1.0           # candidate poses closer than this are the same pose
POSE_TOL_M = 2.0
IMAGE_TOP_N = 3              # imagery is asked to confirm the best neighbour-ranked poses only
SUPPORT_MARGIN = 0.80        # the runner-up pose must explain less than this share of the best
                             # pose's boundary, or the neighbours have not really chosen
PRINTED_REACH_M = 10.0       # a placed parcel the sheet names as a neighbour must be this close


def colour_of(row):
    rms = row.get("boundary_rms_m")
    if row.get("method") == "puvi-only":
        return "red"
    if row.get("printed_far"):
        # the sheet names a neighbour that is on the ground, and the pose does not reach it
        return "red"
    if rms is not None and rms > MAX_RESIDUAL_M:
        # 43A landed 971 m out on 2026-09-18 with a 27.6 m boundary residual and still read amber,
        # because a high image share alone was enough. A placement that disagrees with its own
        # neighbours by metres is wrong however well it matches field edges.
        return "red"
    if row.get("ambiguous") and rms is None:
        return "red"
    allow = SIGMA_AUTO + (row.get("anchor_rms_m") or 0.0)
    anchored = rms is not None and rms <= allow
    strong_image = (row.get("share") or 0.0) >= GREEN_SHARE and row.get("observable")
    weak_image = (row.get("share") or 0.0) >= AMBER_SHARE and row.get("observable")
    # 42A on 2026-09-19: two poses explained 348 m and 347 m of boundary (a slide along the
    # railway strip), the residual was 1.5 m, and the parcel was 6.4 m out. A tie between poses
    # means the neighbours did not choose; imagery may break it, geometry alone may not.
    best, nxt = row.get("support_m") or 0.0, row.get("support_next_m") or 0.0
    decisive = best > 0 and nxt <= SUPPORT_MARGIN * best
    if anchored and (strong_image or (decisive and (row.get("n_neighbours") or 0) >= 2)):
        return "green"
    if anchored or strong_image or weak_image:
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


def neighbour_chains(village, survey, pose, placed, anchor_map, gate_m=CHAIN_GATE_M):
    """Gated shared-boundary observations for one candidate pose.

    `match.common_chains` matches by edge length alone, so a rectangle can match a congruent
    rectangle anywhere in the village. On 2026-09-18 that put 43A's boundary onto chains 78 m away
    and the adjustment carried the parcel 971 m off. A chain is only this parcel's boundary if the
    two rings already nearly coincide there, so any chain whose paired points sit further apart
    than the search window is dropped.

    Returns (observations, supported_length_m, partners).
    """
    theta, t = pose
    v_s = [tuple(q) for q in fit.transform_points(
        sheets.outline(sheets.load_sheet(village, survey)), theta, t)]
    obs, supported, partners = [], 0.0, []
    for other, (th_o, t_o) in placed.items():
        if other == survey or neighbours.are_neighbours(village, survey, other) is False:
            continue
        v_o = [tuple(q) for q in fit.transform_points(
            sheets.outline(sheets.load_sheet(village, other)), th_o, t_o)]
        for chain in match.common_chains(v_s, v_o):
            gap = float(np.median([np.linalg.norm(np.asarray(pa) - np.asarray(pb))
                                   for pa, pb in chain.pairs]))
            if gap > gate_m:
                continue
            sigma = SIGMA_AUTO + (anchor_map[other].rms if other in anchor_map else 0.0)
            # the solver applies each survey's pose to its observation points, so they must be
            # handed over in that survey's own sheet metres. Passing the ground points here moved
            # a correctly chosen 46B 72 m on 2026-09-19.
            local = match.Chain(chain.length, [
                (tuple(fit.unapply_pose(pa, theta, t)[0]), tuple(fit.unapply_pose(pb, th_o, t_o)[0]))
                for pa, pb in chain.pairs], chain.n_pairs, chain.corners)
            obs += match.chain_observations(local, survey, other, sigma)
            supported += chain.length
            partners.append(other)
            break                     # the best surviving chain per neighbour
    return obs, supported, sorted(set(partners), key=paths.survey_sort_key)


def pose_candidates(village, survey, placed, min_len=CHAIN_MIN_M, max_rms=CHAIN_FIT_RMS):
    """Every pose implied by laying this sheet's outline onto a placed neighbour's boundary.

    The chainage walk compares edge lengths and turn angles only, so it can be run between the
    sheet in its own metres and a neighbour already on the ground. Each matching run of edges is
    then a set of point pairs (sheet, ground) whose rigid fit IS the pose. Nothing here uses Puvi:
    the position comes from the team's own placements, which is the rule for this corridor.
    """
    local = sheets.outline(sheets.load_sheet(village, survey))
    out = []
    for other, (th_o, t_o) in placed.items():
        if other == survey or neighbours.are_neighbours(village, survey, other) is False:
            continue
        ground = [tuple(q) for q in fit.transform_points(
            sheets.outline(sheets.load_sheet(village, other)), th_o, t_o)]
        for chain in match.common_chains(local, ground, min_len=min_len):
            P = np.array([pa for pa, _pb in chain.pairs], float)
            Q = np.array([pb for _pa, pb in chain.pairs], float)
            theta, tt, rms, _mx = fit.rigid_fit(P, Q)
            if rms > max_rms:
                continue
            out.append({"theta": float(theta), "t": np.asarray(tt, float), "from": other,
                        "chain_m": chain.length, "fit_rms": float(rms)})
    return out


def dedupe_poses(cands, tol_deg=POSE_TOL_DEG, tol_m=POSE_TOL_M):
    """Collapse candidates that describe the same placement; `votes` counts how many chains agree.

    40B produced 148 candidates on 2026-09-19 (a strip with many equal edges) and the satellite
    search alone then took over nine minutes. Most were the same pose from different chains.
    """
    out = []
    for c in sorted(cands, key=lambda c: -c["chain_m"]):
        for u in out:
            d_th = abs(((c["theta"] - u["theta"]) + 180) % 360 - 180)
            if d_th <= tol_deg and _landing_gap(c, u) <= tol_m:
                u["votes"] += 1
                u["chain_m"] = max(u["chain_m"], c["chain_m"])
                break
        else:
            out.append(dict(c, votes=1))
    return out


def _landing_gap(c, u):
    # t is the translation after rotating about the UTM origin, so two poses that differ by a
    # fraction of a degree differ in t by kilometres; compare where a sheet point actually lands
    probe = np.array([[50.0, 50.0]])
    a = fit.transform_points(probe, c["theta"], c["t"])[0]
    b = fit.transform_points(probe, u["theta"], u["t"])[0]
    return float(np.linalg.norm(a - b))


def placed_bodies(village, placed):
    return {other: unary_union([fit.apply_pose(g, th_o, t_o)
                                for _p, g in sheets.load_sheet(village, other)])
            for other, (th_o, t_o) in placed.items()}


def _overlap(village, survey, pose, bodies):
    """Share of the smaller parcel that this pose steals from something already on the ground."""
    body = unary_union([fit.apply_pose(g, pose[0], pose[1])
                        for _p, g in sheets.load_sheet(village, survey)])
    worst = 0.0
    for other, og in bodies.items():
        if other == survey:
            continue
        inter = body.intersection(og).area
        if inter > 0:
            worst = max(worst, inter / max(min(body.area, og.area), 1.0))
    return worst


def printed_contradictions(village, survey, pose, bodies, reach_m=PRINTED_REACH_M):
    """Placed parcels this sheet names as neighbours that the pose does not come near.

    47A on 2026-09-19 sat 129 m from the team's placement with a 0.75 m residual against one
    long congruent boundary. Its own sheet names 46B, 47B and 48A around it, all placed, and the
    pose touched none of them. The transcription is the cheapest truth we have.
    """
    body = unary_union([fit.apply_pose(g, pose[0], pose[1])
                        for _p, g in sheets.load_sheet(village, survey)])
    far = []
    for other in neighbours.printed(village, survey):
        if other in bodies and other != survey and body.distance(bodies[other]) > reach_m:
            far.append(other)
    return sorted(far, key=paths.survey_sort_key)


def rank_poses(village, survey, candidates, placed, anchor_map, bodies=None):
    """Candidates that do not eat a neighbour, best-supported first: [(support_m, partners, pose)].

    A pose that lands away from a placed parcel the sheet itself names is kept only as a last
    resort, ranked below every pose that honours the printed neighbours.
    """
    bodies = bodies if bodies is not None else placed_bodies(village, placed)
    scored = []
    for pose in candidates:
        if _overlap(village, survey, pose, bodies) > MAX_OVERLAP_SHARE:
            continue
        _obs, supported, partners = neighbour_chains(village, survey, pose, placed, anchor_map)
        far = printed_contradictions(village, survey, pose, bodies)
        scored.append((supported, len(partners), pose, len(far)))
    scored.sort(key=lambda s: (s[3], -s[0], -s[1]))
    return [(s[0], s[1], s[2]) for s in scored], [s[3] for s in scored]


def choose_pose(village, survey, candidates, placed, anchor_map, image=None, ranked=None):
    """Score candidate poses by how much shared boundary they explain without eating a neighbour."""
    scored = list(ranked) if ranked is not None else rank_poses(village, survey, candidates, placed, anchor_map)[0]
    if image is not None and all(p is not image for _s, _n, p in scored):
        scored += rank_poses(village, survey, [image], placed, anchor_map)[0]
        scored.sort(key=lambda s: (-s[0], -s[1]))
    if not scored:
        if image is not None:
            return image, "no pose agrees with the neighbours; pose from imagery"
        return None, "no pose agrees with the neighbours"
    best = scored[0]
    if len(scored) > 1 and scored[1][0] > 0:
        note = "neighbour support %.0f m, next best %.0f m" % (best[0], scored[1][0])
    else:
        note = "neighbour support %.0f m" % best[0]
    return best[2], note


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

    bodies = placed_bodies(village, placed)
    for survey in todo:
        row = {"survey": survey, "run": stamp, "method": "puvi-only", "ambiguous": False,
               "notes": "", "n_neighbours": 0, "file": paths.output_path(village, survey).name}
        cands = dedupe_poses(pose_candidates(village, survey, placed))
        poses = [(c["theta"], c["t"]) for c in cands]
        row["n_candidates"] = len(poses)
        ranked, far_counts = rank_poses(village, survey, poses, placed, anchor_map, bodies)
        image_pose = None
        if index is not None and ranked:
            pts, brg = sheets.sample_outline(sheets.outline(sheets.load_sheet(village, survey)), 1.0)
            res = align.search(pts, brg, index, [p for _s, _n, p in ranked[:IMAGE_TOP_N]])
            row.update({"share": round(res.share, 3), "margin": round(res.margin, 3),
                        "observable": res.observable, "ambiguous": res.ambiguous, "notes": res.note})
            if not res.ambiguous and res.observable:
                image_pose = (res.theta, res.t)
        pose, why = choose_pose(village, survey, poses, placed, anchor_map, image=image_pose,
                                ranked=ranked)
        row["support_m"] = round(ranked[0][0], 1) if ranked else 0.0
        row["support_next_m"] = round(ranked[1][0], 1) if len(ranked) > 1 else 0.0
        if pose is not None:
            far = printed_contradictions(village, survey, pose, bodies)
            row["printed_far"] = ",".join(far)
            if far:
                row["notes"] = " | ".join(x for x in (row.get("notes"),
                                                       "does not reach printed neighbour(s) " + ",".join(far)) if x)
        row["notes"] = " | ".join(x for x in (row.get("notes"), why) if x)
        if pose is not None:
            row["method"] = "image" if (image_pose is not None and pose is image_pose) else "neighbour"
        if pose is None:
            rows.append(dict(row, colour="red", confidence=0, status="waiting",
                             notes=(row["notes"] + " | no placed neighbour shares a boundary").strip(" |")))
            continue
        row["heading_deg"] = round(pose[0], 3)
        placed[survey] = pose
        bodies[survey] = unary_union([fit.apply_pose(g, pose[0], pose[1])
                                      for _p, g in sheets.load_sheet(village, survey)])
        rows.append(row)

    # block adjustment: anchors fixed, everything else free
    free = {s: placed[s] for s in todo if s in placed}
    pair_obs, line_obs, gcp_obs, priors = [], [], [], []
    for s in free:
        priors.append(fit.PosePrior(s, free[s][0], tuple(free[s][1]), SIGMA_START, 10.0))
        obs, _supported, _partners = neighbour_chains(village, s, free[s], placed, anchor_map)
        pair_obs += obs
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
        row["heading_deg"] = round(theta, 3)          # the adjusted heading, not the start guess
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
