"""Joint pose solve for a hand-georeferenced village: shapes from the sheets, positions solved.

    python -m autogeoref.pose_solver 35_04_074 --validate
    python -m autogeoref.pose_solver 35_04_074 --apply

Root cause settled 2026-09-23: the sheets of adjacent surveys agree on their shared printed
boundaries, so the fabric tiles by construction and every overlap or gap is a POSITION error.
Shapes are therefore immutable here: each parcel is regenerated from its sheet at the solved
pose, and nothing is ever clipped or stretched.

The naive solve slid parcels hundreds of metres (51: 542 m) because chains up to
CHAIN_GATE_M = 40 m apart were admitted with sigma 0.30 against a 4 m prior. The guards:

  * GATE_M = 5: a chain is only this parcel's boundary if the rings nearly coincide already.
  * balanced weights: observation sigma 0.75 m against a 2.5 m position prior.
  * freeze-and-resolve: any parcel the solve wants to move further than MOVE_CAP_M is frozen
    at its current pose, its observations dropped, and the solve repeated - a bad chain is cut
    out of the caravan instead of dragging it.
  * leave-one-out validation against the QC-passed parcels, on copies, before --apply will
    touch a single file; --apply refuses when validation fails.
"""
import argparse
import datetime
import logging
import shutil
import sys

import geopandas as gpd
import numpy as np

from . import anchors, engine, fit, fmb_on_base, match, paths, review, sheets, visible

log = logging.getLogger(__name__)

GATE_M = 5.0
OBS_SIGMA = 0.75
PRIOR_POS_M = 2.5
PRIOR_HEAD_DEG = 1.5
MOVE_CAP_M = 6.0
LOO_MEDIAN_M = 1.5     # a kicked parcel must come back at least this close (kick is 3.0 m)
LOO_MAX_M = 3.0        # and never end further out than it started
RAIL = {"582", "583", "584", "607", "608", "609", "610", "611", "612", "613", "614",
        "526B", "569A", "581", "192193207208"}


def village_poses(village):
    """(poses, skipped): the current pose of every manual parcel that has a sheet."""
    poses, skipped = {}, []
    vd = paths.vector_dir(village)
    for f in sorted(vd.glob("*_parcels_modified.gpkg")):
        s = f.name.split("_parcels_modified")[0]
        try:
            sheets.outline(sheets.load_sheet(village, s))
        except Exception:
            skipped.append(s)
            continue
        pose = anchors.pose_from_geometry(village, s, f)
        if pose is None:
            skipped.append(s)
            continue
        theta, t, _rms, _scale, _aniso, _n = pose
        poses[s] = (float(theta), np.asarray(t, float))
    return poses, skipped


NEAR_M = 5.0           # acquisition: a printed neighbour's boundary within this is the seam
SAMPLE_STEP_M = 4.0
CORNER_REACH_M = 3.5
CORNER_SIGMA = 0.6     # corners pin the along-edge direction the boundary samples cannot


def _ground_poly(village, s, pose):
    from shapely.geometry import Polygon
    th, t = pose
    return Polygon(fit.transform_points(sheets.outline(sheets.load_sheet(village, s)), th, t))


def observations(village, placed, free):
    """Boundary-sample and corner correspondences between PRINTED PDF neighbours only.

    Edge-run chains were the wrong instrument here: they match congruent edges anywhere in the
    village (median gaps 24-150 m on 2026-09-23) and miss the true seam whenever the two sheets
    subdivide the shared line differently. Nearest-boundary samples across the seam constrain the
    across-edge direction; paired corners constrain along-edge and heading.
    """
    from shapely.geometry import Point
    from . import neighbours as nbmod
    polys = {s: _ground_poly(village, s, placed[s]) for s in placed}
    obs = []
    for s in free:
        printed = set(nbmod.printed(village, s))
        for other in placed:
            if other == s or (other not in printed and s not in set(nbmod.printed(village, other))):
                continue
            pa, pb = polys[s], polys[other]
            if pa.distance(pb) > NEAR_M:
                continue
            ring = pa.exterior
            n = max(6, int(ring.length / SAMPLE_STEP_M))
            segs = fmb_on_base._segments(pb) if hasattr(fmb_on_base, "_segments") else None
            from .fmb_stage2 import _segments, _nearest_segment
            segs = _segments(pb)
            taken = 0
            for k in range(n):
                pt = np.array(ring.interpolate(k * ring.length / n).coords[0])[:2]
                if pb.exterior.distance(Point(*pt)) > NEAR_M:
                    continue
                (q1, q2), _d = _nearest_segment(pt, segs)
                obs.append(fit.PairLineObs(s, tuple(fit.unapply_pose(pt, *placed[s])[0]),
                                           other, tuple(fit.unapply_pose(q1, *placed[other])[0]),
                                           tuple(fit.unapply_pose(q2, *placed[other])[0]),
                                           OBS_SIGMA))
                taken += 1
                if taken >= 10:
                    break
            ca = fmb_on_base.corners(pa)
            cb = fmb_on_base.corners(pb)
            if len(ca) and len(cb):
                d = np.linalg.norm(ca[:, None, :] - cb[None, :, :], axis=2)
                used = set()
                for i in np.argsort(d.min(axis=1)):
                    j = int(np.argmin(d[i]))
                    if j in used or d[i, j] > CORNER_REACH_M:
                        continue
                    used.add(j)
                    obs.append(match.PairObs(s, other,
                                             tuple(fit.unapply_pose(ca[i], *placed[s])[0]),
                                             tuple(fit.unapply_pose(cb[j], *placed[other])[0]),
                                             CORNER_SIGMA))
    return obs


def solve(village, poses=None, hold_out=(), extra_free=()):
    """Freeze-and-resolve joint solve. Returns (solution poses, frozen surveys)."""
    if poses is None:
        poses, _ = village_poses(village)
    placed = dict(poses)
    free_names = {s for s in placed
                  if (s not in RAIL or s in extra_free) and s not in hold_out} | set(hold_out)
    frozen = set()
    for _round in range(6):
        free = {s: placed[s] for s in free_names - frozen}
        fixed = {s: placed[s] for s in placed if s not in free}
        if not free:
            return {}, frozen
        obs = observations(village, placed, set(free))
        obs = [o for o in obs if o.a not in frozen and o.b not in frozen]
        # production conditions for everyone, the held-out parcel included: its prior sits
        # at its (wrong) current pose with the same sigma - the test is whether the seams
        # move it toward the truth, not whether it can be found with no information at all
        priors = [fit.PosePrior(s, th, tuple(t), PRIOR_POS_M, PRIOR_HEAD_DEG)
                  for s, (th, t) in free.items()]
        point_obs = [o for o in obs if isinstance(o, match.PairObs)]
        line_obs = [o for o in obs if isinstance(o, fit.PairLineObs)]
        sol = fit.block_adjust(free, fixed, point_obs, [], [], priors, pair_line_obs=line_obs)
        over = {s for s, r in sol.items() if r["shift_m"] > MOVE_CAP_M and s not in hold_out}
        if not over:
            return sol, frozen
        frozen |= over
        log.info("froze %s: the solve asked for more than %.1f m", sorted(over), MOVE_CAP_M)
    return sol, frozen


PERTURB_M = 3.0
PERTURB_DEG = 2.0


def validate(village):
    """Hold each trusted parcel out, PERTURB it, and demand the chains pull it back.

    Without the perturbation the test is circular: a settled fabric already satisfies its own
    chains and every hold-out scores a meaningless 0.00.
    """
    poses, _ = village_poses(village)
    approved = [r["survey"] for r in review.read_status(village)
                if r.get("status") == "approved" and r["survey"] in poses]
    if not approved:
        approved = [s for s in poses if s not in RAIL]      # a fully hand-placed village
    rows = []
    for k, s in enumerate(approved):
        bump = np.array([PERTURB_M, 0.0]) if k % 2 == 0 else np.array([0.0, -PERTURB_M])
        pert = dict(poses)
        pert[s] = (poses[s][0] + (PERTURB_DEG if k % 2 == 0 else -PERTURB_DEG),
                   poses[s][1] + bump)
        sol, frozen = solve(village, poses=pert, hold_out=(s,))
        n_obs = len([o for o in observations(village, pert, {s})])
        if s not in sol:
            rows.append({"survey": s, "err_m": float("nan"), "note": "frozen or unsolved"})
            continue
        err = float(np.linalg.norm(np.asarray(sol[s]["t"], float) - poses[s][1]))   # vs TRUE pose
        head = abs(((sol[s]["theta"] - poses[s][0]) + 180) % 360 - 180)
        rows.append({"survey": s, "err_m": round(err, 2), "head_deg": round(head, 2),
                     "note": "%d obs" % n_obs})
    errs = [r["err_m"] for r in rows if r["err_m"] == r["err_m"]]
    verdict = {"parcels": len(rows), "solved": len(errs),
               "median_m": round(float(np.median(errs)), 2) if errs else None,
               "max_m": round(float(max(errs)), 2) if errs else None,
               "adopt": bool(errs) and float(np.median(errs)) <= LOO_MEDIAN_M
                        and float(max(errs)) <= LOO_MAX_M}
    return rows, verdict


def apply(village):
    """Disabled: it rewrote Akash's hand files in place (2026-09-24 recovery).

    A pose change now goes to the tool's own _auto folder and reaches his files only through
    an approved install. The solver itself stays available through --validate.
    """
    raise RuntimeError("pose_solver --apply is disabled: it would rewrite hand files; "
                       "write proposals to paths.auto_dir() and install them after approval")


def _apply_retired(village):
    _rows, verdict = validate(village)
    log.info("validation: %s", verdict)
    if not verdict["adopt"]:
        log.error("validation FAILED - nothing applied")
        return None
    poses, _ = village_poses(village)
    sol, frozen = solve(village, poses=poses)
    vd = paths.vector_dir(village)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    bdir = paths.logs_dir() / ("backup_%s_posesolve_%s" % (village, ts))
    bdir.mkdir(parents=True, exist_ok=True)
    moved = {}
    for s, r in sol.items():
        f = vd / ("%s_parcels_modified.gpkg" % s)
        visible.refuse_tool_write(f)
        shutil.copy2(f, bdir / f.name)
        plots = sheets.load_sheet(village, s)
        rows = [{"plot_no": str((p or {}).get("plot_no", "")), "survey_no": s,
                 "geometry": fit.apply_pose(g, r["theta"], r["t"])} for p, g in plots]
        gpd.GeoDataFrame(rows, geometry="geometry", crs=32644).to_file(
            f, layer="parcels", driver="GPKG")
        moved[s] = round(r["shift_m"], 2)
    log.info("applied %d parcels (backup %s); frozen: %s", len(moved), bdir.name, sorted(frozen))
    return {"moved": moved, "frozen": sorted(frozen), "backup": str(bdir)}


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="pose_solver", description=__doc__.splitlines()[0])
    ap.add_argument("village")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args(argv)
    if args.validate or not args.apply:
        rows, verdict = validate(args.village)
        for r in rows:
            print("  %-5s err %.2f m %s" % (r["survey"], r["err_m"], r["note"]))
        print("verdict:", verdict)
        return 0 if verdict["adopt"] else 2
    out = apply(args.village)
    return 0 if out else 2


if __name__ == "__main__":
    sys.exit(main())
