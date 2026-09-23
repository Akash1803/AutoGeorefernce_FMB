"""Place each FMB sheet on the fine-tuned Puvi base.

    python -m autogeoref.fmb_on_base --corridor --base "<Puvi_Vector_Finetunned.geojson>"

The sheet gives the geometry, which is the legal record and carries the subdivision plots; the
base gives the position. Each sheet is fitted onto its base parcel by rotation and shift only,
with scale fixed at 1, so every printed length survives.

Akash's rule, 2026-09-22: try the shape match first; where a parcel's shape would have to change to
fit, keep the base's shape instead and record the disagreement. Control points are taken from the
base layer's own corners, paired with the sheet's, and written in the QGIS Georeferencer format so
a placement can be dragged by hand.

Measured on the three villages his placements anchor: shape match 0.86 to 0.91, area ratio 1.00,
and the placed parcel inherits the base's error exactly (2.7 to 4.2 m). This method cannot be more
accurate than the base, and it is not less.
"""
import argparse
import csv
import datetime
import logging
import math
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.validation import make_valid

from . import anchors, paths, shift_puvi

log = logging.getLogger(__name__)

UTM = 32644
MATCH_MIN = 0.60          # below this the sheet and the base describe different ground
AREA_BAND = (0.80, 1.25)  # and so does an area ratio outside this
CORNER_TURN_DEG = 20.0    # a vertex that turns less than this is not a corner anyone can point at
CORNER_MIN_ARM_M = 1.0
GCP_REACH_M = 5.0         # a base corner and a sheet corner this close are the same corner
AMBIGUOUS_TURN = 0.05     # two fits this close in score are a parcel that fits two ways


def _valid(geom):
    if geom is None or geom.is_empty:
        return None
    return geom if geom.is_valid else make_valid(geom)


def sheet_of(village, survey):
    """The sheet's plots and its outline, in sheet metres. Returns (GeoDataFrame, outline)."""
    f = paths.vector_dir(village) / ("%s_parcels.geojson" % survey)
    if not f.exists():
        return None, None
    try:
        g = gpd.read_file(f)
    except Exception as exc:
        log.warning("%s %s: sheet unreadable (%s)", village, survey, exc)
        return None, None
    geoms = [_valid(x) for x in g.geometry]
    good = [x for x in geoms if x is not None]
    if not good:
        return None, None
    return g, unary_union(good)


def corners(geom, turn_deg=CORNER_TURN_DEG, arm_m=CORNER_MIN_ARM_M):
    """Vertices where the outline really turns: the points a person could put a pin in."""
    out = []
    for poly in shift_puvi._polygons_of(geom):
        # drop the closing point: with it, vertex 0's previous neighbour is itself, which made the
        # first corner of every parcel invisible
        ring = np.array(poly.exterior.coords)[:-1, :2]
        n = len(ring)
        for i in range(n):
            a, b, c = ring[i - 1], ring[i], ring[(i + 1) % n]
            v1, v2 = a - b, c - b
            l1, l2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
            if l1 < arm_m or l2 < arm_m:
                continue
            ang = math.degrees(math.acos(max(-1.0, min(1.0, float(v1 @ v2) / (l1 * l2)))))
            if abs(180.0 - ang) >= turn_deg:
                out.append(b)
    return np.array(out) if out else np.zeros((0, 2))


def pair_corners(placed_geom, base_geom, reach=GCP_REACH_M):
    """Corners of the placed sheet matched to the base's. Returns [(sheet_xy, base_xy, dist)]."""
    cs, cb = corners(placed_geom), corners(base_geom)
    if not len(cs) or not len(cb):
        return []
    d = np.linalg.norm(cs[:, None, :] - cb[None, :, :], axis=2)
    taken, out = set(), []
    for i in np.argsort(d.min(axis=1)):
        j = int(np.argmin(d[i]))
        if j in taken or d[i, j] > reach:
            continue
        taken.add(j)
        out.append((cs[i], cb[j], float(d[i, j])))
    return out


def place(sheet_outline, base_geom):
    """Rotation and shift putting the sheet over its base parcel. Returns a dict, or None."""
    if sheet_outline is None or base_geom is None or sheet_outline.is_empty or base_geom.is_empty:
        return None
    ang, shift, score = shift_puvi.rigid_fit(sheet_outline, base_geom)
    flip = (ang + 180.0) % 360.0
    moved = affinity.translate(sheet_outline, shift[0], shift[1])
    origin = Point(*base_geom.centroid.coords[0])
    other = affinity.rotate(moved, flip, origin=origin)
    inter = other.intersection(base_geom).area
    flip_score = inter / (other.area + base_geom.area - inter) if other.area else 0.0
    return {"angle": float(ang), "shift": np.asarray(shift, float), "match": float(score),
            "area_ratio": float(sheet_outline.area / base_geom.area) if base_geom.area else float("nan"),
            "ambiguous": bool(abs(flip_score - score) <= AMBIGUOUS_TURN)}


def apply_pose(geom, pose, base_geom):
    """Put one plot of the sheet where the pose says."""
    if geom is None:
        return None
    moved = affinity.translate(geom, pose["shift"][0], pose["shift"][1])
    return _valid(affinity.rotate(moved, pose["angle"], origin=Point(*base_geom.centroid.coords[0])))


REPORT_COLUMNS = ["village_code", "village", "survey_no", "geometry_source", "shape_match",
                  "area_ratio", "boundary_gap_m", "plots_expected", "plots_kept",
                  "self_overlap_sqm", "gap_to_base_sqm", "gcp_points", "gcp_residual_m",
                  "ambiguous_turn", "evidence", "flag", "run_at"]


def conservation(placed_plots, sheet_count, base_geom, placed_outline):
    """What the placement kept and what it lost, in counts and square metres.

    Taken from the tool Akash shared: a placement that quietly drops plots, overlaps itself or
    leaves the base's ground uncovered should say so, not read as clean.
    """
    kept = [g for g in placed_plots if g is not None and not g.is_empty]
    overlap = 0.0
    for i in range(len(kept)):
        for j in range(i + 1, len(kept)):
            inter = kept[i].intersection(kept[j])
            if not inter.is_empty:
                overlap += inter.area
    body = unary_union(kept) if kept else None
    gap = float(base_geom.difference(body).area) if body is not None else float(base_geom.area)
    snap = float(placed_outline.hausdorff_distance(base_geom.boundary)) if placed_outline is not None else float("nan")
    return {"plots_expected": int(sheet_count), "plots_kept": len(kept),
            "self_overlap_sqm": round(overlap, 1), "gap_to_base_sqm": round(gap, 1),
            "boundary_gap_m": round(snap, 2)}


def run_village(village, base, out_dir, hand=None, stamp=None):
    """Place every sheet of one village. Returns (plot rows, report rows, unplaced rows)."""
    hand = hand or {}
    stamp = stamp or datetime.datetime.now().isoformat(timespec="seconds")
    name = str(base["village"].iloc[0]) if "village" in base.columns and len(base) else village
    plots, report, unplaced = [], [], []
    gcp_dir = out_dir / "gcp"

    for _, row in base.iterrows():
        survey = str(row["survey_no"])
        base_geom = _valid(row.geometry)
        common = {"village_code": village, "village": name, "survey_no": survey,
                  "evidence": row.get("evidence", ""), "run_at": stamp}
        if base_geom is None:
            unplaced.append({**common, "reason": "the base parcel has no geometry"})
            continue

        if survey in hand and hand[survey] is not None:
            plots.append({**common, "geometry_source": "your placement", "plot_no": "",
                          "shape_match": "", "area_ratio": "", "gcp_points": "",
                          "geometry": hand[survey]})
            report.append({**common, "geometry_source": "your placement", "shape_match": "",
                           "area_ratio": "", "gcp_points": "", "gcp_residual_m": "",
                           "plots_expected": 1, "plots_kept": 1, "ambiguous_turn": "", "flag": ""})
            continue

        sheet_gdf, outline = sheet_of(village, survey)
        if outline is None:
            unplaced.append({**common, "reason": "no converted sheet"})
            plots.append({**common, "geometry_source": "puvi base", "plot_no": "",
                          "shape_match": "", "area_ratio": "", "gcp_points": "",
                          "geometry": base_geom})
            report.append({**common, "geometry_source": "puvi base", "shape_match": "",
                           "area_ratio": "", "gcp_points": "", "gcp_residual_m": "",
                           "plots_expected": 1, "plots_kept": 1, "ambiguous_turn": "",
                           "flag": "no sheet"})
            continue

        pose = place(outline, base_geom)
        if pose is None:
            unplaced.append({**common, "reason": "the sheet could not be fitted"})
            continue

        agrees = pose["match"] >= MATCH_MIN and AREA_BAND[0] <= pose["area_ratio"] <= AREA_BAND[1]
        flag = ""
        if not agrees:
            flag = ("shape match %.2f" % pose["match"] if pose["match"] < MATCH_MIN
                    else "area ratio %.2f" % pose["area_ratio"])
        if pose["ambiguous"]:
            flag = (flag + "; fits two ways").strip("; ")

        if agrees:
            placed_outline = apply_pose(outline, pose, base_geom)
            pairs = pair_corners(placed_outline, base_geom)
            for _i, srow in sheet_gdf.iterrows():
                g = apply_pose(_valid(srow.geometry), pose, base_geom)
                if g is None:
                    continue
                plots.append({**common, "geometry_source": "fmb sheet",
                              "plot_no": str(srow.get("plot_no", "")),
                              "shape_match": round(pose["match"], 3),
                              "area_ratio": round(pose["area_ratio"], 3),
                              "gcp_points": len(pairs), "geometry": g})
            cons = conservation([p["geometry"] for p in plots[-len(sheet_gdf):]] if len(sheet_gdf) else [],
                                len(sheet_gdf), base_geom, placed_outline)
            if cons["plots_kept"] < cons["plots_expected"]:
                flag = (flag + "; %d of %d plots kept" % (cons["plots_kept"], cons["plots_expected"])).strip("; ")
            if pairs:
                gcp_dir.mkdir(parents=True, exist_ok=True)
                rows = []
                for sheet_xy, base_xy, _d in pairs:
                    # the sheet coordinate is where the point sits on the drawing, before placing
                    back = affinity.rotate(Point(*sheet_xy), -pose["angle"],
                                           origin=Point(*base_geom.centroid.coords[0]))
                    back = affinity.translate(back, -pose["shift"][0], -pose["shift"][1])
                    rows.append((float(base_xy[0]), float(base_xy[1]), float(back.x), float(-back.y)))
                anchors.write_points(gcp_dir / ("%s_%s_parcels.geojson.points" % (village, survey)),
                                     rows, "EPSG:%d" % UTM)
            report.append({**common, "geometry_source": "fmb sheet",
                           "shape_match": round(pose["match"], 3),
                           "area_ratio": round(pose["area_ratio"], 3),
                           "gcp_points": len(pairs),
                           "gcp_residual_m": round(float(np.median([p[2] for p in pairs])), 2) if pairs else "",
                           "ambiguous_turn": int(pose["ambiguous"]), "flag": flag, **cons})
        else:
            # the shapes disagree: the base's shape stands, as agreed
            plots.append({**common, "geometry_source": "puvi base", "plot_no": "",
                          "shape_match": round(pose["match"], 3),
                          "area_ratio": round(pose["area_ratio"], 3), "gcp_points": 0,
                          "geometry": base_geom})
            report.append({**common, "geometry_source": "puvi base",
                           "shape_match": round(pose["match"], 3),
                           "area_ratio": round(pose["area_ratio"], 3), "gcp_points": 0,
                           "gcp_residual_m": "", "plots_expected": len(sheet_gdf), "plots_kept": 1,
                           "ambiguous_turn": int(pose["ambiguous"]), "flag": flag})
    return plots, report, unplaced


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="fmb_on_base", description=__doc__.splitlines()[0])
    ap.add_argument("--base", required=True, help="the fine-tuned Puvi layer")
    ap.add_argument("--corridor", action="store_true", help="every village in the base")
    ap.add_argument("--villages", nargs="*", default=[], help="village codes, if not the whole base")
    ap.add_argument("--control", action="append", default=[],
                    help="extra file holding your own placements (repeatable)")
    ap.add_argument("--out", default="", help="output folder")
    args = ap.parse_args(argv)

    base = gpd.read_file(args.base).to_crs(UTM)
    villages = args.villages or (sorted(set(base["village_code"].astype(str))) if args.corridor else [])
    if not villages:
        ap.error("name some villages, or pass --corridor")
    stamp = datetime.datetime.now().isoformat(timespec="seconds")
    out_dir = Path(args.out) if args.out else (
        paths.PROJECT / ("FMB_on_Puvi_%s" % datetime.date.today().strftime("%Y%m%d")))
    out_dir.mkdir(parents=True, exist_ok=True)

    all_plots, all_report, all_unplaced = [], [], []
    for v in villages:
        sub = base[base["village_code"].astype(str) == v]
        if not len(sub):
            continue
        hand = {k: _valid(g) for k, g in shift_puvi.control_parcels(v, args.control).items()}
        plots, report, unplaced = run_village(v, sub, out_dir, hand=hand, stamp=stamp)
        all_plots += plots
        all_report += report
        all_unplaced += unplaced
        log.info("%s: %d plot(s) from %d survey(s); %d on the base's shape, %d unplaced",
                 v, len(plots), len(report),
                 sum(1 for r in report if r["geometry_source"] == "puvi base"), len(unplaced))

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
    if all_unplaced:
        pd.DataFrame(all_unplaced).drop(columns=["geometry"], errors="ignore").to_csv(
            out_dir / "unplaced.csv", index=False)

    src = pd.Series([r["geometry_source"] for r in all_report]).value_counts().to_dict()
    print("\n%d plots from %d surveys over %d villages" % (len(g), len(all_report), len(villages)))
    print("geometry: %s" % src)
    print("flagged for your eye: %d" % sum(1 for r in all_report if r.get("flag")))
    print("written to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
