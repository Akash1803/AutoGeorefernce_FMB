"""Make neighbouring placed sheets agree: stage 2 of the FMB-on-base SOP.

    python -m autogeoref.fmb_stage2 --gate --base "<Puvi_Vector_Finetunned.geojson>"
    python -m autogeoref.fmb_stage2 --corridor --base "<Puvi_Vector_Finetunned.geojson>"

Stage 1 places each sheet on its own base parcel, so neighbouring sheets meet only as well as the
base does. Here each placed sheet is adjusted against the sheets it touches: who is a neighbour
comes from the Puvi fabric, the boundary geometry comes from the sheets, and the adjustment is a
rotation and shift per parcel, scale still locked at 1. Akash's placements, the parcels keeping
the base's shape, and the railway land never move.

The gate, agreed 2026-09-22: this stage ships only if, measured on the parcels Akash placed by
hand, it moves parcels closer to his placements. Closing gaps is not the test; being right is.
The gate run holds none of his parcels fixed, so nothing it is scored on ever leaks in - which
also makes it conservative: the real run has his 53 as extra anchors.
"""
import argparse
import csv
import datetime
import logging
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import Point
from shapely.ops import unary_union

from . import engine, fit, fmb_on_base, paths, shift_puvi, shift_score
from .fit import PairLineObs, PosePrior
from .match import LineObs

log = logging.getLogger(__name__)

UTM = 32644
TOUCH_M = 1.0          # base parcels this close are neighbours in the fabric
SAMPLE_STEP_M = 5.0    # boundary samples along a placed sheet's outline
NEAR_M = 8.0           # a sample further than this from the neighbour is a disagreement, not a seam
EDGE_SIGMA_M = 1.0
PRIOR_POS_M = 3.0      # the base's own accuracy: how far a parcel may stray from stage 1
PRIOR_HEAD_DEG = 1.5
MAX_OBS_PER_PAIR = 12
MAX_MOVE_M = 6.0       # a solve asking for more than this is closing a disagreement; refused


# --------------------------------------------------------------------------- geometry helpers

def _boundary_samples(geom, step=SAMPLE_STEP_M):
    out = []
    for poly in shift_puvi._polygons_of(geom):
        ring = poly.exterior
        n = max(4, int(ring.length / step))
        out += [np.array(ring.interpolate(i * ring.length / n).coords[0])[:2] for i in range(n)]
    return out


def _segments(geom):
    segs = []
    for poly in shift_puvi._polygons_of(geom):
        c = np.array(poly.exterior.coords)[:, :2]
        segs += [(c[i], c[i + 1]) for i in range(len(c) - 1)]
    return segs


def _nearest_segment(p, segs):
    best, best_d = None, None
    for q1, q2 in segs:
        d = q2 - q1
        L2 = float(d @ d)
        t = max(0.0, min(1.0, float((p - q1) @ d) / L2)) if L2 else 0.0
        dist = float(np.linalg.norm(q1 + t * d - p))
        if best_d is None or dist < best_d:
            best, best_d = (q1, q2), dist
    return best, best_d


def neighbour_pairs(base, touch=TOUCH_M):
    """Unordered survey pairs whose BASE parcels touch: the Puvi fabric names the neighbours."""
    rows = [(str(r["survey_no"]), fmb_on_base._valid(r.geometry)) for _, r in base.iterrows()]
    rows = [(k, g) for k, g in rows if g is not None]
    out = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if rows[i][0] != rows[j][0] and rows[i][1].distance(rows[j][1]) <= touch:
                out.append((rows[i][0], rows[j][0]))
    return sorted(set(out))


# --------------------------------------------------------------------------- the adjustment

def _observations(geoms, centres, pairs, free):
    """Boundary samples of each free sheet tied to the neighbour's nearest edge."""
    pair_line, line = [], []
    for a, b in pairs:
        for s, t in ((a, b), (b, a)):
            if s not in free or s not in geoms or t not in geoms:
                continue
            t_boundary = geoms[t].boundary
            segs = _segments(geoms[t])
            near = [p for p in _boundary_samples(geoms[s])
                    if t_boundary.distance(Point(*p)) <= NEAR_M]
            if len(near) > MAX_OBS_PER_PAIR:
                near = [near[i] for i in
                        np.linspace(0, len(near) - 1, MAX_OBS_PER_PAIR).astype(int)]
            for p in near:
                (q1, q2), dist = _nearest_segment(p, segs)
                if t in free:
                    pair_line.append(PairLineObs(s, tuple(p - centres[s]), t,
                                                 tuple(q1 - centres[t]), tuple(q2 - centres[t]),
                                                 EDGE_SIGMA_M))
                else:
                    line.append(LineObs(s, tuple(p - centres[s]), (tuple(q1), tuple(q2)),
                                        EDGE_SIGMA_M, dist))
    return pair_line, line


def adjust(geoms, pairs, free):
    """A small rotation + shift per free survey so neighbouring sheets meet.

    `geoms` are the stage-1 placed outlines in UTM; `free` names the surveys allowed to move.
    Returns {survey: delta} for the free surveys; everything else is an anchor.
    """
    free = {s for s in free if s in geoms}
    centres = {s: np.array(g.centroid.coords[0]) for s, g in geoms.items()}
    if not free:
        return {}
    pair_line, line = _observations(geoms, centres, pairs, free)
    free_pose = {s: (0.0, tuple(centres[s])) for s in free}
    fixed_pose = {s: (0.0, tuple(centres[s])) for s in geoms if s not in free}
    priors = [PosePrior(s, 0.0, tuple(centres[s]), PRIOR_POS_M, PRIOR_HEAD_DEG) for s in free]
    sol = fit.block_adjust(free_pose, fixed_pose, [], line, [], priors, pair_line_obs=pair_line)

    out = {}
    for s in free:
        c = centres[s]
        theta = sol[s]["theta"] if s in sol else 0.0
        t = np.asarray(sol[s]["t"], float) if s in sol else c
        delta = {"theta": float(theta), "t": t, "origin": c,
                 "moved_m": float(np.linalg.norm(t - c)), "flag": ""}
        # the guard: how far the worst vertex would travel
        worst = 0.0
        for p in np.array(unary_union(shift_puvi._polygons_of(geoms[s])).envelope.exterior.coords)[:, :2]:
            worst = max(worst, float(np.linalg.norm(
                fit.transform_points([p - c], theta, t)[0] - p)))
        if worst > MAX_MOVE_M:
            delta = {"theta": 0.0, "t": c, "origin": c, "moved_m": 0.0,
                     "flag": "refused: solve asked for %.1f m" % worst}
        out[s] = delta
    return out


def apply_delta(geom, delta):
    if geom is None:
        return None
    c, t = delta["origin"], delta["t"]
    g = affinity.rotate(geom, delta["theta"], origin=(float(c[0]), float(c[1])))
    return affinity.translate(g, float(t[0] - c[0]), float(t[1] - c[1]))


# --------------------------------------------------------------------------- the gate

def gate_verdict(rows, tol=0.05):
    """Adopt only if the parcels Akash placed end up closer to his placements."""
    if not rows:
        return {"adopt": False, "parcels": 0, "reason": "nothing could be measured"}
    b = np.array([r["before_m"] for r in rows])
    a = np.array([r["after_m"] for r in rows])
    rmse_b, rmse_a = float(np.sqrt(np.mean(b ** 2))), float(np.sqrt(np.mean(a ** 2)))
    med_b, med_a = float(np.median(b)), float(np.median(a))
    adopt = (rmse_a < rmse_b + tol and med_a < med_b + tol
             and (rmse_a < rmse_b or med_a < med_b))
    return {"adopt": bool(adopt), "parcels": len(rows),
            "rmse_before_m": round(rmse_b, 2), "rmse_after_m": round(rmse_a, 2),
            "median_before_m": round(med_b, 2), "median_after_m": round(med_a, 2),
            "better": int((a < b - 0.2).sum()), "worse": int((a > b + 0.2).sum())}


def _stage1_outlines(village, base_sub):
    """Each survey's stage-1 placed outline and where its geometry came from, hand control unused."""
    geoms, sources = {}, {}
    for _, row in base_sub.iterrows():
        survey = str(row["survey_no"])
        base_geom = fmb_on_base._valid(row.geometry)
        if base_geom is None:
            continue
        _gdf, outline = fmb_on_base.sheet_of(village, survey)
        if outline is None:
            geoms[survey], sources[survey] = base_geom, "no sheet"
            continue
        pose = fmb_on_base.place(outline, base_geom)
        agrees = (pose is not None and pose["match"] >= fmb_on_base.MATCH_MIN
                  and fmb_on_base.AREA_BAND[0] <= pose["area_ratio"] <= fmb_on_base.AREA_BAND[1])
        if agrees:
            geoms[survey] = fmb_on_base.apply_pose(outline, pose, base_geom)
            sources[survey] = "fmb sheet"
        else:
            geoms[survey], sources[survey] = base_geom, "puvi base"
    return geoms, sources


def _free_of(village, sources):
    rail = engine._rail_parcels(village)
    return {s for s, src in sources.items() if src == "fmb sheet" and s not in rail}


def gate(base, control_files=()):
    """Score stage 2 against Akash's placements, none of them fixed, so nothing leaks."""
    rows = []
    for v in shift_score.SCORED_VILLAGES:
        truth = shift_score.truth_for(v, control_files if v == "35_04_052" else ())
        if not truth:
            continue
        sub = base[base["village_code"].astype(str) == v]
        geoms, sources = _stage1_outlines(v, sub)
        deltas = adjust(geoms, neighbour_pairs(sub), _free_of(v, sources))
        for k, (t, _s) in truth.items():
            if sources.get(k) != "fmb sheet":
                continue
            c = np.array(geoms[k].centroid.coords[0])
            after = np.asarray(deltas[k]["t"], float) if k in deltas else c
            rows.append({"village": v, "survey": k,
                         "before_m": float(np.linalg.norm(t - c)),
                         "after_m": float(np.linalg.norm(t - after)),
                         "moved_m": deltas.get(k, {}).get("moved_m", 0.0),
                         "flag": deltas.get(k, {}).get("flag", "")})
    return rows, gate_verdict(rows)


# --------------------------------------------------------------------------- the full run

def run_village(village, base_sub, out_dir, hand, stamp):
    """Stage 1 then the neighbour adjustment for one village. Returns (plots, report, unplaced, moves)."""
    plots, report, unplaced = fmb_on_base.run_village(village, base_sub, out_dir,
                                                     hand=hand, stamp=stamp)
    by_survey = {}
    for p in plots:
        by_survey.setdefault(str(p["survey_no"]), []).append(p)
    geoms = {s: fmb_on_base._valid(unary_union([p["geometry"] for p in rows]))
             for s, rows in by_survey.items()}
    geoms = {s: g for s, g in geoms.items() if g is not None}
    sources = {str(r["survey_no"]): r["geometry_source"] for r in report}
    free = _free_of(village, sources)
    deltas = adjust(geoms, neighbour_pairs(base_sub), free)

    moves = []
    for s, d in sorted(deltas.items()):
        for p in by_survey.get(s, []):
            p["geometry"] = apply_delta(p["geometry"], d)
            p["moved_m"] = round(d["moved_m"], 2)
        moves.append({"village_code": village, "survey_no": s, "moved_m": round(d["moved_m"], 2),
                      "dtheta_deg": round(d["theta"], 3), "flag": d["flag"]})
    for p in plots:
        p.setdefault("moved_m", 0.0)
    for r in report:
        d = deltas.get(str(r["survey_no"]))
        r["moved_m"] = round(d["moved_m"], 2) if d else 0.0
        if d and d["flag"]:
            r["flag"] = (str(r.get("flag", "")) + "; " + d["flag"]).strip("; ")
    return plots, report, unplaced, moves


REPORT_COLUMNS = fmb_on_base.REPORT_COLUMNS + ["moved_m"]


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="fmb_stage2", description=__doc__.splitlines()[0])
    ap.add_argument("--base", required=True, help="the fine-tuned Puvi layer")
    ap.add_argument("--corridor", action="store_true", help="every village in the base")
    ap.add_argument("--villages", nargs="*", default=[], help="village codes, if not the whole base")
    ap.add_argument("--control", action="append", default=[],
                    help="extra file holding Akash's own placements (repeatable)")
    ap.add_argument("--gate", action="store_true",
                    help="only score the stage against his placements; write nothing else")
    ap.add_argument("--out", default="", help="output folder")
    args = ap.parse_args(argv)

    base = gpd.read_file(args.base).to_crs(UTM)
    out_dir = Path(args.out) if args.out else (
        paths.PROJECT / ("FMB_on_Puvi_%s" % datetime.date.today().strftime("%Y%m%d")))
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.gate:
        rows, verdict = gate(base, tuple(args.control))
        pd.DataFrame(rows).to_csv(out_dir / "stage2_gate.csv", index=False)
        print("\ngate on %(parcels)d of Akash's parcels: RMSE %(rmse_before_m).2f -> "
              "%(rmse_after_m).2f m, median %(median_before_m).2f -> %(median_after_m).2f m, "
              "%(better)d better / %(worse)d worse" % verdict)
        print("ADOPT stage 2" if verdict["adopt"] else "DO NOT ADOPT: stage 1 ships alone")
        return 0 if verdict["adopt"] else 2

    villages = args.villages or (sorted(set(base["village_code"].astype(str))) if args.corridor else [])
    if not villages:
        ap.error("name some villages, or pass --corridor")
    stamp = datetime.datetime.now().isoformat(timespec="seconds")

    all_plots, all_report, all_unplaced, all_moves = [], [], [], []
    for v in villages:
        sub = base[base["village_code"].astype(str) == v]
        if not len(sub):
            continue
        hand = {k: fmb_on_base._valid(g)
                for k, g in shift_puvi.control_parcels(v, args.control).items()}
        plots, report, unplaced, moves = run_village(v, sub, out_dir, hand, stamp)
        all_plots += plots
        all_report += report
        all_unplaced += unplaced
        all_moves += moves
        moved = [m for m in moves if m["moved_m"] > 0.05]
        log.info("%s: %d survey(s), %d adjusted (median %.2f m), %d refused",
                 v, len(report), len(moved),
                 float(np.median([m["moved_m"] for m in moved])) if moved else 0.0,
                 sum(1 for m in moves if m["flag"]))

    if not all_plots:
        log.error("nothing placed")
        return 1
    g = gpd.GeoDataFrame(all_plots, geometry="geometry", crs=UTM)
    g["area_sqm"] = g.geometry.area.round(1)
    g.to_crs(4326).to_file(out_dir / "FMB_parcels_georeferenced.geojson",
                           driver="GeoJSON", COORDINATE_PRECISION=8)
    by_survey = g.dissolve(by=["village_code", "survey_no"], aggfunc="first").reset_index()
    by_survey.to_crs(4326).to_file(out_dir / "FMB_parcels_by_survey.geojson",
                                   driver="GeoJSON", COORDINATE_PRECISION=8)
    with (out_dir / "placement_report.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS)
        w.writeheader()
        for r in all_report:
            w.writerow({c: r.get(c, "") for c in REPORT_COLUMNS})
    pd.DataFrame(all_moves).to_csv(out_dir / "stage2_moves.csv", index=False)
    if all_unplaced:
        pd.DataFrame(all_unplaced).drop(columns=["geometry"], errors="ignore").to_csv(
            out_dir / "unplaced.csv", index=False)

    moved = [m for m in all_moves if m["moved_m"] > 0.05]
    print("\n%d plots from %d surveys; %d surveys adjusted, median move %.2f m, largest %.2f m"
          % (len(g), len(all_report), len(moved),
             float(np.median([m["moved_m"] for m in moved])) if moved else 0.0,
             max((m["moved_m"] for m in moved), default=0.0)))
    print("refused (a disagreement, not a seam): %d" % sum(1 for m in all_moves if m["flag"]))
    print("written to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
