"""The per-village pipeline of spec section 5.

Order: freeze the team's work, pick a start pose for everything else, let the satellite imagery
correct that pose where the evidence is unambiguous, tie the free parcels to their placed
neighbours in one block adjustment, write the files, then report. Anchors are never moved and
never rewritten; Puvi only ever suggests a starting guess and a reporting distance.
"""
import datetime
import math
from pathlib import Path

import geopandas as gpd
import numpy as np
from pyproj import Proj
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from . import (align, anchors, edges as edgemod, files, fit, gcp, match, neighbours,
               paths, raster, review, sheets, topology, visible, window)

ACRE = 0.000247105381
SIGMA_AUTO = 0.30
SIGMA_IMAGE = 1.0
SIGMA_START = 10.0
IMAGE_MAX_SHIFT = 1.5        # imagery may move an anchored pose by at most this
IMAGE_MAX_TURN = 2.0
# Calibrated on the Kizhikaranai leave-one-out of 2026-09-19 (15 parcels, table in spec section 3).
# Satellite edges were observable on only 3 of 15 parcels there, and the one parcel they alone
# would have passed (48A, share 0.81) was 4.1 m out, so the image thresholds are deliberately high.
GREEN_SHARE = 0.85
AMBER_SHARE = 0.45
GREEN_MIN_NEIGHBOURS = 3     # every parcel with three supporting neighbours landed within 2.7 m;
                             # with two, 47B/48A/48B were 4-9 m out and looked no different
DISTINCT_M = 5.0             # a runner-up pose must land this far from the best, or turn
DISTINCT_DEG = 3.0           # this much, to count as a different answer (else it is the same pose)
CHAIN_GATE_M = 40.0          # a chain whose paired points are further apart than the search window
                             # is a congruent boundary somewhere else, not this parcel's boundary
MAX_RESIDUAL_M = 5.0         # above this the placement contradicts its own neighbours
CHAIN_MIN_M = 12.0           # a shared run shorter than this pins nothing
CHAIN_FIT_RMS = 2.0          # metres: how well a chain must fit before its pose is a candidate
ADJUST_GATE_M = 8.0          # chains admitted to the adjustment: correct ones sit within the
                             # anchors' disagreement (about 3 m); wrong ones were 13 m and more
MAX_OVERLAP_SHARE = 0.02     # a pose may not sit on top of a parcel already on the ground
OVERLAP_TOL_M = 3.0          # the team's own placements disagree by up to 3 m; an overlap strip
                             # that thin along a shared edge is not "sitting on top of" anything
POSE_TOL_DEG = 1.0           # candidate poses closer than this are the same pose
POSE_TOL_M = 2.0
IMAGE_TOP_N = 3              # imagery is asked to confirm the best neighbour-ranked poses only
SUPPORT_MARGIN = 0.80        # the runner-up pose must explain less than this share of the best
                             # pose's boundary, or the neighbours have not really chosen
PRINTED_REACH_M = 10.0       # a placed parcel the sheet names as a neighbour must be this close
SIDE_COS = 0.38              # a neighbour must lie within 67.5 degrees of the side the sheet prints it on
LINE_REACH_M = 8.0           # a boundary sample may seek a printed neighbour's edge this far away
LINE_STEP_M = 2.0            # metres between boundary samples that become line observations
LINE_DEG = 15.0              # the neighbour's edge must run within this of the sample's own edge
SIGMA_LINE = 0.60            # metres; many samples per edge, so each one is weighted loosely
ADJUST_ROUNDS = 3            # rebuild the observations from the adjusted poses and solve again
WINDOW_MARGIN_M = 40.0       # satellite window: the placed parcels' extent plus this on every side


def colour_of(row):
    rms = row.get("boundary_rms_m")
    if row.get("method") == "puvi-only":
        return "red"
    far = [x for x in str(row.get("printed_far") or "").split(",") if x]
    contra = row.get("contradictions")
    if contra is None:
        contra = len(far) + (row.get("side_bad") or 0)
    if contra >= 2:
        # the sheets name neighbours and sides that this pose does not honour
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
    n = row.get("n_neighbours") or 0
    if n <= 1 and not (strong_image or weak_image):
        # 40B, 2026-09-19: one neighbour, residual 0.24 m, 16 m out. One chain pins a line, not
        # a parcel, and nothing else vouched for it.
        return "red"
    if contra:
        return "amber" if (anchored or strong_image or weak_image) else "red"
    tied_to_team = (row.get("anchor_partners") if row.get("anchor_partners") is not None else n) >= 1
    if anchored and tied_to_team and ((strong_image and n >= 2) or (decisive and n >= GREEN_MIN_NEIGHBOURS)):
        # 47B read green at 126 m on 2026-09-19 with three "neighbours" that were all placed in the
        # same pass from the same wrong guess; green needs a direct tie to the team's own work
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

    visible.refuse_tool_write(target)
    if target.exists():
        if not files.safe_write(target, writer, allow_release=True):
            raise PermissionError("cannot replace %s (open in QGIS?)" % target)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        writer(target)
    return target


def boundary_line_observations(village, survey, pose, placed, anchor_map, reach_m=LINE_REACH_M):
    """Tie this sheet's boundary samples to the edges of the printed neighbours already placed.

    For each sample (every LINE_STEP_M along the outline, in sheet metres) the nearest edge of a
    printed neighbour's sheet is found in ground coordinates; it must run within LINE_DEG of the
    sample's own edge and lie within reach_m. The observation is stored in both sheets' own
    metres, so it holds whether the neighbour is fixed or still free.
    """
    theta, t = pose
    local_pts, brg = sheets.sample_outline(sheets.outline(sheets.load_sheet(village, survey)), LINE_STEP_M)
    ground = fit.transform_points(local_pts, theta, t)
    out = []
    for other, (th_o, t_o) in placed.items():
        if other == survey or neighbours.are_neighbours(village, survey, other) is not True:
            continue
        v_loc = sheets.outline(sheets.load_sheet(village, other))
        v_gnd = fit.transform_points(v_loc, th_o, t_o)
        sigma = SIGMA_LINE + (anchor_map[other].rms if other in anchor_map else 0.0)
        n = len(v_loc)
        edges_g = [(v_gnd[i], v_gnd[(i + 1) % n]) for i in range(n)]
        for k, p in enumerate(ground):
            sb = (brg[k] + math.radians(theta)) % math.pi
            best = None
            for i, (a, b) in enumerate(edges_g):
                d = b - a
                L = float(np.hypot(d[0], d[1]))
                if L < 1.0:
                    continue
                eb = math.atan2(d[1], d[0]) % math.pi
                diff = abs((eb - sb + math.pi / 2) % math.pi - math.pi / 2)
                if diff > math.radians(LINE_DEG):
                    continue
                # perpendicular distance, but only within the edge's own extent
                s = float(((p - a) @ d) / (L * L))
                if s < -0.05 or s > 1.05:
                    continue
                perp = abs(((p[0] - a[0]) * d[1] - (p[1] - a[1]) * d[0]) / L)
                if best is None or perp < best[0]:
                    best = (perp, i)
            if best is None or best[0] > reach_m:
                continue
            i = best[1]
            out.append(fit.PairLineObs(survey, (float(local_pts[k][0]), float(local_pts[k][1])), other,
                                       (float(v_loc[i][0]), float(v_loc[i][1])),
                                       (float(v_loc[(i + 1) % n][0]), float(v_loc[(i + 1) % n][1])), sigma))
    return out


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


def _overlap(village, survey, pose, bodies, tol_m=OVERLAP_TOL_M):
    """Share of the smaller parcel that this pose steals from something already on the ground.

    Both parcels are eroded by `tol_m` first. On 2026-09-19 the unfiltered test rejected the
    correct pose of 47A (0.1 m from the team's) because a 2 m strip along its 120 m edge with
    171 was 10 % of its area; every thin parcel and all three railway strips lost their correct
    candidate the same way. A wrong pose that really sits on a neighbour survives the erosion.
    """
    body = unary_union([fit.apply_pose(g, pose[0], pose[1])
                        for _p, g in sheets.load_sheet(village, survey)])
    core = body.buffer(-tol_m)
    if core.is_empty:
        core = body                                # a parcel thinner than 2 * tol_m: compare as is
    worst = 0.0
    for other, og in bodies.items():
        if other == survey:
            continue
        oc = og.buffer(-tol_m)
        if oc.is_empty:
            oc = og
        inter = core.intersection(oc).area
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


def side_checks(village, survey, pose, placed, bodies):
    """(honoured, contradicted): the compass sides the sheets print, read both ways.

    The sheet under placement says on which side each printed neighbour lies, and every placed
    neighbour's own sheet says on which side this survey lies. Sides are in the sheet's drawn
    frame, so they turn with the pose. From a single fixed parcel this is the only evidence that
    tells which edge of the anchor a neighbour attaches to: on 2026-09-19, seeded with 48A alone,
    47B took a wrong edge and every later parcel followed it hundreds of metres away.
    """
    body = unary_union([fit.apply_pose(g, pose[0], pose[1])
                        for _p, g in sheets.load_sheet(village, survey)])
    c0 = np.array(body.centroid.coords[0])
    ok = bad = 0

    def rot(theta):
        th = math.radians(theta)
        return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])

    R = rot(pose[0])
    for o in neighbours.printed(village, survey):
        letters = neighbours.side(village, survey, o)
        if o == survey or o not in bodies or not letters:
            continue
        act = np.array(bodies[o].centroid.coords[0]) - c0
        n = float(np.linalg.norm(act))
        if n < 1.0:
            continue
        if float((R @ np.array(neighbours.side_vector(letters))) @ act) / n >= SIDE_COS:
            ok += 1
        else:
            bad += 1
    for o, (th_o, _t_o) in placed.items():
        letters = neighbours.side(village, o, survey)
        if o == survey or o not in bodies or not letters:
            continue
        act = c0 - np.array(bodies[o].centroid.coords[0])
        n = float(np.linalg.norm(act))
        if n < 1.0:
            continue
        if float((rot(th_o) @ np.array(neighbours.side_vector(letters))) @ act) / n >= SIDE_COS:
            ok += 1
        else:
            bad += 1
    return ok, bad


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
        s_ok, s_bad = side_checks(village, survey, pose, placed, bodies)
        scored.append((supported, len(partners), pose, len(far) + s_bad, s_ok))
    # contradictions first (fewest), then support, then how many printed sides agree
    scored.sort(key=lambda s: (s[3], -s[0], -s[4], -s[1]))
    return [(s[0], s[1], s[2]) for s in scored], [s[3] for s in scored]


def _runner_up(village, survey, ranked):
    """Support of the best pose that is a different answer from the winner.

    Many candidates are the same placement reached through different chains and differ by
    centimetres; on 2026-09-19 that made support_next equal support on 13 of 15 parcels and
    nothing could ever be green.
    """
    best = ranked[0][2]
    probe = np.array([[50.0, 50.0]])
    b = fit.transform_points(probe, best[0], best[1])[0]
    for support, _n, pose in ranked[1:]:
        d_th = abs(((pose[0] - best[0]) + 180) % 360 - 180)
        gap = float(np.linalg.norm(fit.transform_points(probe, pose[0], pose[1])[0] - b))
        if gap >= DISTINCT_M or d_th >= DISTINCT_DEG:
            return support
    return 0.0


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


def _place_one(village, survey, stamp, placed, bodies, anchor_map, index, pass_no):
    """Propose a pose for one sheet from what is on the ground now: (row, pose or None)."""
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
    row["support_next_m"] = round(_runner_up(village, survey, ranked), 1) if ranked else 0.0
    if pose is not None:
        far = printed_contradictions(village, survey, pose, bodies)
        s_ok, s_bad = side_checks(village, survey, pose, placed, bodies)
        row["printed_far"] = ",".join(far)
        row["side_ok"], row["side_bad"] = s_ok, s_bad
        row["contradictions"] = len(far) + s_bad
        if far:
            row["notes"] = " | ".join(x for x in (row.get("notes"),
                                                   "does not reach printed neighbour(s) " + ",".join(far)) if x)
        if s_bad:
            row["notes"] = " | ".join(x for x in (row.get("notes"),
                                                   "%d printed side(s) contradicted" % s_bad) if x)
        # Akash, 2026-09-23: the sheet's printed neighbours are checked BEFORE a placement is
        # finalised. A pose the sheet itself argues against is refused, never delivered:
        # 47B went 90 m up the rail strip with 9 contradictions and still reached review.
        if row["contradictions"] >= 2 and row["contradictions"] > (s_ok or 0):
            row["notes"] = " | ".join(x for x in (row.get("notes"),
                                                   "refused: the printed neighbours disagree with every pose") if x)
            return dict(row, colour="red", confidence=0,
                        status="refused: printed neighbours disagree"), None
    row["notes"] = " | ".join(x for x in (row.get("notes"), why) if x)
    if pose is not None:
        row["method"] = "image" if (image_pose is not None and pose is image_pose) else "neighbour"
    if pose is None:
        return dict(row, colour="red", confidence=0, status="waiting"), None
    row["heading_deg"] = round(pose[0], 3)
    row["pass"] = pass_no
    return row, pose


def _commit(village, survey, pose, placed, bodies):
    placed[survey] = pose
    bodies[survey] = unary_union([fit.apply_pose(g, pose[0], pose[1])
                                  for _p, g in sheets.load_sheet(village, survey)])


def _decisive(row):
    best, nxt = row.get("support_m") or 0.0, row.get("support_next_m") or 0.0
    return best > 0 and nxt <= SUPPORT_MARGIN * best


def _place_all(village, todo, stamp, placed, bodies, anchor_map, index):
    """The placement passes: propose every pending sheet, commit the certain ones, refine, repeat.

    `placed` and `bodies` are updated in place. Returns {"rows": [...], "bodies": {survey: body}}
    for the sheets placed here (anchors excluded).
    """
    rows = []
    # Placement grows outward in passes: a sheet whose neighbours are not on the ground yet is
    # retried after they are. With one fixed parcel (48A on 2026-09-19) the whole village is
    # reached in a handful of passes; a sheet nothing ever reaches ends as "waiting".
    # Within a pass every pending sheet is proposed against the same ground, then only the
    # certain ones are committed: no contradiction with the printed neighbours and sides, and a
    # clear winner over the runner-up. If none is certain, the single best-supported proposal is
    # committed so the village keeps growing; the colour rule will mark it for review.
    pending, pass_no = list(todo), 0
    while pending:
        pass_no += 1
        proposals = []
        for survey in pending:
            row, pose = _place_one(village, survey, stamp, placed, bodies, anchor_map, index, pass_no)
            if pose is not None:
                proposals.append((row, pose))
        certain = [(r, p) for r, p in proposals if (r.get("contradictions") or 0) == 0 and _decisive(r)]
        chosen = certain or sorted(proposals, key=lambda rp: ((rp[0].get("contradictions") or 0),
                                                              -(rp[0].get("support_m") or 0)))[:1]
        chosen.sort(key=lambda rp: -(rp[0].get("support_m") or 0))
        placed_this_pass = []
        for row, pose in chosen:
            if _overlap(village, row["survey"], pose, bodies) > MAX_OVERLAP_SHARE:
                continue                       # collides with a sheet committed earlier this pass
            _commit(village, row["survey"], pose, placed, bodies)
            row["certain"] = (row, pose) in certain
            rows.append(row)
            placed_this_pass.append(row["survey"])
        # refinement: a sheet committed early in a pass never saw the sheets committed after it.
        # 47B (2026-09-20) was proposed against 48A alone, took the one short shared edge and
        # came out 9.5 degrees off; its long edge with 171, committed in the same pass, would
        # have pinned it. Re-propose each newly placed sheet against everything now on the
        # ground and keep the new pose when more boundary supports it.
        for survey in placed_this_pass:
            others = {k: p for k, p in placed.items() if k != survey}
            other_bodies = {k: b for k, b in bodies.items() if k != survey}
            row2, pose2 = _place_one(village, survey, stamp, others, other_bodies, anchor_map, index, pass_no)
            old_row = next(r for r in rows if r["survey"] == survey)
            if pose2 is not None and (row2.get("contradictions") or 0) <= (old_row.get("contradictions") or 0)                     and (row2.get("support_m") or 0) > (old_row.get("support_m") or 0) + 1.0:
                row2["certain"] = old_row.get("certain")
                row2["notes"] = " | ".join(x for x in (row2.get("notes"), "refined after pass %d" % pass_no) if x)
                rows[rows.index(old_row)] = row2
                _commit(village, survey, pose2, placed, bodies)
        pending = [s for s in pending if s not in placed_this_pass]
        if not placed_this_pass:
            for survey in pending:
                rows.append({"survey": survey, "run": stamp, "method": "puvi-only", "ambiguous": False,
                             "n_neighbours": 0, "file": paths.output_path(village, survey).name,
                             "colour": "red", "confidence": 0, "status": "waiting",
                             "notes": "no placed neighbour shares a boundary (after %d passes)" % pass_no})
            break
    return {"rows": rows, "bodies": {s: bodies[s] for s in bodies if s in todo}}

def run(village, do_raster=False, do_topology=False, do_review=True, project=None, only=None,
        raster_bounds=None, use_imagery=True):
    """Steps 1-13 of spec section 5. Returns a summary dict; details go to georef_status.csv.

    `only`: place just these surveys (the ring around a seed); everything else is left alone.
    `raster_bounds`: export the satellite window for these EPSG:32644 bounds instead of the
    whole village. A small window keeps the tile reprojection error local (Akash, 2026-09-20).
    """
    run_id = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    tool_fps = review.tool_written_fingerprints(village)
    anchor_map = anchors.load_anchors(village, tool_fps)
    conflicts = anchors.anchor_conflicts(village, anchor_map)
    superseded = {s: paths.manual_files(village, s)[1:] for s in anchor_map}
    anchors.write_anchors_csv(village, anchor_map, superseded, conflicts)

    todo = [s for s in paths.surveys_with_sheets(village) if s not in anchor_map]
    if only is not None:
        only = {str(s) for s in only}
        todo = [s for s in todo if s in only]
    placed = {s: (a.theta, a.t) for s, a in anchor_map.items()}
    rows = []

    auto_window = raster_bounds == "auto"
    if auto_window:
        # Sized from where the ring lands, not from the seed: around 48A alone the window missed
        # the far end of the 300 m strip 171 and its GCPs there had no imagery (Akash, 2026-09-20).
        pre = _place_all(village, todo, stamp, dict(placed), placed_bodies(village, placed), anchor_map, None)
        geoms = [b for s, b in pre["bodies"].items()] + [anchors.placed_geometry(village, a) for a in anchor_map.values()]
        geoms = [g for g in geoms if g is not None and not g.is_empty]
        b = unary_union(geoms).bounds
        raster_bounds = (b[0] - WINDOW_MARGIN_M, b[1] - WINDOW_MARGIN_M, b[2] + WINDOW_MARGIN_M, b[3] + WINDOW_MARGIN_M)
        do_raster = True
    if not use_imagery:
        # neighbours only: no window is fetched and no edge index is built. Thirukatchur's 19
        # anchors span the whole village and the automatic window for them ran past ten minutes.
        do_raster, raster_bounds = False, None
    if todo and use_imagery and (do_raster or raster.pinned(village) is None):
        if raster_bounds is not None:
            raster.export(village, tuple(raster_bounds), run_id=run_id)
        else:
            geoms = [anchors.placed_geometry(village, a) for a in anchor_map.values()]
            geoms = [g for g in geoms if g is not None and not g.is_empty]
            if geoms:
                b = unary_union(geoms).bounds
                raster.export(village, (b[0] - 100, b[1] - 100, b[2] + 100, b[3] + 100), run_id=run_id)

    index = None
    if todo and use_imagery and raster.pinned(village) is not None:
        index = edgemod.SegmentIndex(edgemod.detect(raster.pinned(village)["path"], min_len_m=3.0))

    bodies = placed_bodies(village, placed)
    result = _place_all(village, todo, stamp, placed, bodies, anchor_map, index)
    rows.extend(result["rows"])
    if index is not None:
        # every placed parcel must lie inside the satellite window, or its GCPs are guesses
        pin = raster.pinned(village)
        rb = pin.get("bounds_32644") if pin else None
        if rb:
            for row in rows:
                b = bodies.get(row["survey"])
                if b is not None and not (b.bounds[0] >= rb[0] and b.bounds[1] >= rb[1]
                                          and b.bounds[2] <= rb[2] and b.bounds[3] <= rb[3]):
                    row["notes"] = " | ".join(x for x in (row.get("notes"), "outside the satellite window") if x)
                    row["in_window"] = False
                else:
                    row["in_window"] = True

    # block adjustment: anchors fixed, everything else free. Observations are rebuilt from the
    # adjusted poses and the solve repeated, so an edge that was out of reach at first (47B's far
    # end 13 m from 171's line) is caught once the first round has turned the parcel.
    free = {s: placed[s] for s in todo if s in placed}
    fixed_poses = {s: (a.theta, a.t) for s, a in anchor_map.items()}
    line_obs, gcp_obs = [], []
    adjusted, pair_obs = {}, []
    for _round in range(ADJUST_ROUNDS if free else 0):
        current = dict(fixed_poses, **{s: (p[0], p[1]) for s, p in free.items()})
        pair_obs, pair_line_obs, priors = [], [], []
        for s in free:
            priors.append(fit.PosePrior(s, free[s][0], tuple(free[s][1]), SIGMA_START, 10.0))
            obs, _supported, _partners = neighbour_chains(village, s, free[s], current, anchor_map,
                                                          gate_m=ADJUST_GATE_M)
            pair_obs += obs
            pair_line_obs += boundary_line_observations(village, s, free[s], current, anchor_map)
        adjusted = fit.block_adjust(free, fixed_poses, pair_obs, line_obs, gcp_obs, priors,
                                    pair_line_obs=pair_line_obs)
        free = {s: (adjusted[s]["theta"], adjusted[s]["t"]) for s in adjusted}
    for s in adjusted:
        adjusted[s]["line_obs"] = sum(1 for o in pair_line_obs if o.a == s)
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
        row["pose_tx"], row["pose_ty"] = round(float(t[0]), 4), round(float(t[1]), 4)
        mine = [v for k, v in residuals.items() if s in k]
        row["boundary_rms_m"] = round(min(v["rms"] for v in mine), 3) if mine else None
        partners = sorted({k[0] if k[1] == s else k[1] for k in residuals if s in k},
                          key=paths.survey_sort_key)
        row["n_neighbours"] = len(partners)
        row["neighbours"] = ",".join(partners)
        row["anchor_rms_m"] = max([anchor_map[n].rms for n in partners if n in anchor_map] or [0.0])
        row["anchor_partners"] = sum(1 for n in partners if n in anchor_map)
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
        # the team's parcels take part as the team drew them (affine, vertex-edited), not as the
        # rigid sheet at the recovered pose: on 2026-09-20 47B was clipped clean against the
        # rigid 48A and still overlapped Akash's real 48A by 19.5 m2
        geoms.update({s: team_geometry(village, a) for s, a in anchor_map.items()})
        rail = {s for s in geoms if s in _rail_parcels(village)}
        order = sorted(written, key=lambda s: (next((int(r.get("pass") or 99) for r in rows if r["survey"] == s), 99),
                                               paths.survey_sort_key(s)))
        edited = write_topology(village, written, adjusted, {s: geoms[s] for s in anchor_map}, order, rail)
        for row in rows:
            if row["survey"] in edited:
                row["notes"] = " | ".join(x for x in (row.get("notes"), "topology: " + edited[row["survey"]]) if x)
                row["fp_placed"] = anchors.fingerprint(paths.output_path(village, row["survey"]))
                if "topology cut" in edited[row["survey"]] and row.get("colour") != "red":
                    # a parcel the resolver had to carve was mis-placed; the carving is not a fix
                    row["colour"], row["confidence"] = "red", 0

    for row in rows:
        if row.get("fp_placed"):
            row["fp_final"] = row["fp_placed"]
    review.write_status(village, rows)
    if do_review:
        review.build_group(village, rows, project)
    return {"village": village, "anchors": len(anchor_map), "placed": len(written),
            "waiting": sum(1 for r in rows if r.get("status") == "waiting"),
            "conflicts": len(conflicts), "run": run_id, "rows": rows}


def team_geometry(village, anchor):
    """The hand-placed parcel exactly as the team saved it, dissolved to one body."""
    try:
        hand = visible.read_hand(anchor.file)
        if hand.crs is not None and hand.crs.to_epsg() != 32644:
            hand = hand.to_crs(32644)
        body = unary_union([g.buffer(0) for g in hand.geometry if g is not None and not g.is_empty])
        if not body.is_empty:
            return body
    except Exception:
        pass
    return anchors.placed_geometry(village, anchor)


def write_topology(village, written, adjusted, anchor_bodies, order, rail=()):
    """Rewrite the parcels layer of each tool-written file with its topology-resolved plots.

    Only the parcels layer changes; the edges layer stays the rigid sheet, because the printed
    FMB lengths are the point of it. Anchors are never in `written`, so the team's files are never
    touched. Returns {survey: summary} for the surveys whose geometry changed.
    """
    rigid = {}
    for s in written:
        theta, t = adjusted[s]["theta"], adjusted[s]["t"]
        rigid[s] = [(props, fit.apply_pose(g, theta, t)) for props, g in sheets.load_sheet(village, s)]
    resolved, report = topology.resolve(rigid, anchor_bodies, order, rail=rail, anchors=set(anchor_bodies))
    topology.report_rail_conflicts(village, report)
    edited = {}
    for s in written:
        new_parts = resolved.get(s) or [(p, g, "") for p, g in rigid[s]]
        new_parts, guard_notes = topology.sanity(rigid[s], new_parts)
        if any("came apart" in n or "topology cut" in n for n in guard_notes):
            # plot-level conform tore the survey, or clipping carved it: try it as one body against
            # what is settled (small moves only; a real mis-placement still overlaps and is refused)
            settled_now = dict(anchor_bodies)
            settled_now.update({o: unary_union([g for _p, g, _n in resolved[o]]) for o in order
                                if o != s and o in resolved and order.index(o) < order.index(s)})
            alt, why = topology.conform_as_body(s, rigid[s], settled_now, anchors=set(anchor_bodies),
                                                anchor_tol=topology.ANCHOR_CONFORM_TOL * 1.25)
            if alt is not None:
                alt, alt_notes = topology.sanity(rigid[s], alt)
                if not any("rigid-clipped" in n or "came apart" in n for n in alt_notes):
                    new_parts, guard_notes = alt, [why] + alt_notes
        notes = sorted({n for _p, _g, n in new_parts if n}) + guard_notes
        if s in report.get("conformed", {}):
            notes.append("conform max %.2f m" % report["conformed"][s]["max_move_m"])
        if not notes:
            continue
        target = paths.output_path(village, s)
        visible.refuse_tool_write(target)
        old = gpd.read_file(target, layer="parcels")

        def pkey(pid, pno):
            # None, NaN and '' are the same missing plot number; ids compare as text
            missing = pno is None or pno == "" or (isinstance(pno, float) and pno != pno)
            return (str(pid), "" if missing else str(pno))

        geoms_by_id = {pkey(props.get("poly_id"), props.get("plot_no")): (g, n) for props, g, n in new_parts}
        new_geom, topo, unmatched = [], [], 0
        for i, (_, r) in enumerate(old.iterrows()):
            hit = geoms_by_id.get(pkey(r.get("poly_id"), r.get("plot_no")))
            if hit is None and len(old) == len(new_parts):
                hit = (new_parts[i][1], new_parts[i][2])      # same sheet order as write_parcels
            if hit is None:
                unmatched += 1
                hit = (r.geometry, "unmatched")
            new_geom.append(hit[0]); topo.append(hit[1] or "clean")
        if unmatched:
            notes.append("%d plot(s) unmatched" % unmatched)
        old = old.set_geometry(new_geom)
        old["topo_edit"] = topo
        old["area_sqm"] = old.geometry.area.round(3)
        old["perimeter_m"] = old.geometry.length.round(3)

        def writer(path, frame=old, src=target):
            frame.to_file(path, layer="parcels", driver="GPKG")
            gpd.read_file(src, layer="edges").to_file(path, layer="edges", driver="GPKG")

        if files.safe_write(target, writer, allow_release=True):
            edited[s] = ",".join(notes)
    return edited


def pose_from_edges(village, survey, gpkg_path):
    """Exact pose of a tool-written file, from its edges layer.

    The parcels layer may have been topology-edited, and a rigid fit to it was 0.5 degrees off
    for 47B on 2026-09-21 (a metre over its length). The edges layer is the rigid sheet as
    placed, so polygonising each plot's edges and fitting to that returns the pose to the mm.
    """
    from shapely.ops import polygonize
    ed = gpd.read_file(gpkg_path, layer="edges")
    polys = []
    for (pid, pn), grp in ed.groupby(["poly_id", "plot_no"], dropna=False):
        pg = list(polygonize(unary_union(list(grp.geometry))))
        if pg:
            polys.append({"poly_id": pid, "plot_no": pn, "geometry": max(pg, key=lambda q: q.area)})
    tmp = Path(gpkg_path).with_name("_pose_from_edges_%s.gpkg" % survey)
    gpd.GeoDataFrame(polys, geometry="geometry", crs="EPSG:32644").to_file(tmp, driver="GPKG")
    try:
        pose = anchors.pose_from_geometry(village, survey, tmp)
    finally:
        tmp.unlink(missing_ok=True)
    return {"theta": float(pose[0]), "t": np.asarray(pose[1], float)}


def retopology(village, surveys=None):
    """Re-run the topology stage on files the tool already wrote, without placing anything again.

    The pose comes from the status row (pose_tx/pose_ty written since 2026-09-21) or, for older
    files, from a rigid fit of the sheet to the written geometry. Anchors take part as the team's
    geometry and are never written. Returns {survey: summary}.
    """
    rows = review.read_status(village)
    tool_fps = review.tool_written_fingerprints(village)
    anchor_map = anchors.load_anchors(village, tool_fps)
    placed_rows = [r for r in rows if r.get("status") == "placed" and (surveys is None or r["survey"] in surveys)]
    adjusted = {}
    for r in placed_rows:
        s = r["survey"]
        if r.get("pose_tx") and r.get("pose_ty"):
            adjusted[s] = {"theta": float(r["heading_deg"]), "t": np.array([float(r["pose_tx"]), float(r["pose_ty"])])}
        else:
            adjusted[s] = pose_from_edges(village, s, paths.output_path(village, s))
    written = list(adjusted)
    geoms = {s: unary_union([fit.apply_pose(g, adjusted[s]["theta"], adjusted[s]["t"])
                             for _p, g in sheets.load_sheet(village, s)]) for s in written}
    geoms.update({s: team_geometry(village, a) for s, a in anchor_map.items()})
    rail = {s for s in geoms if s in _rail_parcels(village)}
    order = sorted(written, key=lambda s: (next((int(r.get("pass") or 99) for r in placed_rows if r["survey"] == s), 99),
                                           paths.survey_sort_key(s)))
    edited = write_topology(village, written, adjusted, {s: geoms[s] for s in anchor_map}, order, rail)
    for r in rows:
        if r["survey"] in edited:
            base = (r.get("notes") or "").split(" | topology:")[0]
            r["notes"] = base + " | topology: " + edited[r["survey"]]
            r["fp_placed"] = r["fp_final"] = anchors.fingerprint(paths.output_path(village, r["survey"]))
    review.write_status(village, rows)
    return edited


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
