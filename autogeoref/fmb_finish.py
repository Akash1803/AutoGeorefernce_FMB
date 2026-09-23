"""Finish the FMB placement: the topology pass and the QC fields (stages 3 and 4 of the SOP).

    python -m autogeoref.fmb_finish --corridor --base "<Puvi_Vector_Finetunned.geojson>"

Stage 1 places each sheet on its own base parcel; stage 2 (neighbour agreement) failed its gate on
2026-09-23 and is not applied. What remains is the corridor's own topology pass: overlaps between
placed surveys are clipped out of the yielding side, thin seams are filled into the parcel that
borders them most, and no edit may empty a parcel or take more than a quarter of it. Akash's
placements are never edited, and railway land keeps its shape with conflicts reported.

Every plot then carries its evidence: geometry_source, shape_match, area_ratio, gcp_points from
stage 1, and topo_note saying exactly what this pass did to it ("" for the untouched majority,
which keep the sheet's drawing to the millimetre).
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
from shapely.ops import unary_union

from . import engine, fmb_on_base, paths, shift_puvi, shift_score, topology

log = logging.getLogger(__name__)

UTM = 32644
GAP_TOL_M = 1.20        # the corridor's own seam width: wider is a disagreement, not a sliver
MAX_FILL_M2 = 25.0      # the SOP's cap: a bigger "gap" is land the buffer does not hold


def finish_village(village, plots, report, rail=None):
    """Run the topology pass over one village's placed plots. Returns (plots, report, topo)."""
    by_survey = {}
    for p in plots:
        by_survey.setdefault(str(p["survey_no"]), []).append(p)
    bodies, sources = {}, {}
    for s, rows in by_survey.items():
        g = fmb_on_base._valid(unary_union([r["geometry"] for r in rows]))
        if g is not None and not g.is_empty:
            bodies[s] = g
            sources[s] = rows[0].get("geometry_source", "")
    movable = {s for s in bodies if sources.get(s) != "your placement"}
    rail = engine._rail_parcels(village) if rail is None else rail

    fixed, topo = topology.fix(bodies, movable, rail, gap_tol=GAP_TOL_M,
                               max_fill_m2=MAX_FILL_M2)
    edited = {r[0] for r in topo["clips"]} | {r[2] for r in topo["fills"]}

    for s in sorted(edited):
        rows = by_survey.get(s, [])
        if not rows:
            continue
        parts = [(r, r["geometry"]) for r in rows]
        for (r, _g), (_r, new_g, note) in zip(parts,
                                              topology.apply_to_parts(parts, fixed[s], bodies[s])):
            r["geometry"] = new_g
            r["topo_note"] = note
    for p in plots:
        p.setdefault("topo_note", "")

    note_of = {}
    for key, _what, other, area in topo["clips"]:
        note_of[key] = (note_of.get(key, "") + "; clipped %.1f m2 against %s" % (area, other)).strip("; ")
    for area, _into, key, kind in topo["fills"]:
        note_of[key] = (note_of.get(key, "") + "; filled %.1f m2 (%s)" % (area, kind)).strip("; ")
    for key, why, taken, reason in topo["refused"]:
        note_of[key] = (note_of.get(key, "") + "; refused %s (%s)" % (why, reason)).strip("; ")
    for r in report:
        r["topo_note"] = note_of.get(str(r["survey_no"]), "")
    return plots, report, topo


def cross_village_overlap(plots):
    """How much land two villages' placed sheets both claim, in m2 per village pair."""
    bodies = {}
    for p in plots:
        bodies.setdefault(str(p["village_code"]), []).append(p["geometry"])
    bodies = {v: fmb_on_base._valid(unary_union(g)) for v, g in bodies.items()}
    out = {}
    keys = sorted(k for k, g in bodies.items() if g is not None)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            inter = bodies[a].intersection(bodies[b])
            if not inter.is_empty and inter.area > 1.0:
                out[(a, b)] = round(float(inter.area), 1)
    return out


def _write_safe(geom):
    """Valid in the written CRS, not only in metres: snap, then repair what remains."""
    g = shift_puvi._snap_or_keep(geom)
    if g is not None and not g.is_empty and not g.is_valid:
        g = g.buffer(0)
    return g


REPORT_COLUMNS = fmb_on_base.REPORT_COLUMNS + ["topo_note"]


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="fmb_finish", description=__doc__.splitlines()[0])
    ap.add_argument("--base", required=True, help="the fine-tuned Puvi layer")
    ap.add_argument("--corridor", action="store_true", help="every village in the base")
    ap.add_argument("--villages", nargs="*", default=[], help="village codes, if not the whole base")
    ap.add_argument("--control", action="append", default=[],
                    help="extra file holding Akash's own placements (repeatable)")
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
    topo_rows = []
    for v in villages:
        sub = base[base["village_code"].astype(str) == v]
        if not len(sub):
            continue
        hand = {k: fmb_on_base._valid(g)
                for k, g in shift_puvi.control_parcels(v, args.control).items()}
        plots, report, unplaced = fmb_on_base.run_village(v, sub, out_dir, hand=hand, stamp=stamp)
        plots, report, topo = finish_village(v, plots, report)
        all_plots += plots
        all_report += report
        all_unplaced += unplaced
        for key, _what, other, area in topo["clips"]:
            topo_rows.append({"village_code": v, "survey_no": key, "action": "clipped",
                              "against": other, "sqm": area})
        for area, _into, key, kind in topo["fills"]:
            topo_rows.append({"village_code": v, "survey_no": key, "action": "filled",
                              "against": kind, "sqm": area})
        for key, why, taken, reason in topo["refused"]:
            topo_rows.append({"village_code": v, "survey_no": key, "action": "refused",
                              "against": "%s: %s" % (why, reason), "sqm": round(taken, 1)})
        for a, b, area in topo["rail_conflicts"]:
            topo_rows.append({"village_code": v, "survey_no": a, "action": "rail conflict",
                              "against": b, "sqm": area})
        log.info("%s: %d clip(s), %d fill(s), %d refused, %d rail conflict(s)",
                 v, len(topo["clips"]), len(topo["fills"]),
                 len(topo["refused"]), len(topo["rail_conflicts"]))

    if not all_plots:
        log.error("nothing placed")
        return 1
    across = cross_village_overlap(all_plots)
    for (a, b), area in sorted(across.items()):
        log.info("village overlap %s / %s: %.1f m2 (left as drawn; the base owns the boundary)", a, b, area)

    g = gpd.GeoDataFrame(all_plots, geometry="geometry", crs=UTM)
    g["area_sqm"] = g.geometry.area.round(1)
    g = g.to_crs(4326)
    g["geometry"] = [_write_safe(x) for x in g.geometry]
    g.to_file(out_dir / "FMB_parcels_georeferenced.geojson",
              driver="GeoJSON", COORDINATE_PRECISION=8)
    by_survey = g.dissolve(by=["village_code", "survey_no"], aggfunc="first").reset_index()
    by_survey["geometry"] = [_write_safe(x) for x in by_survey.geometry]
    by_survey.to_file(out_dir / "FMB_parcels_by_survey.geojson",
                      driver="GeoJSON", COORDINATE_PRECISION=8)
    with (out_dir / "placement_report.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS)
        w.writeheader()
        for r in all_report:
            w.writerow({c: r.get(c, "") for c in REPORT_COLUMNS})
    pd.DataFrame(topo_rows).to_csv(out_dir / "topology_report.csv", index=False)
    if all_unplaced:
        pd.DataFrame(all_unplaced).drop(columns=["geometry"], errors="ignore").to_csv(
            out_dir / "unplaced.csv", index=False)

    src = pd.Series([r["geometry_source"] for r in all_report]).value_counts().to_dict()
    touched = sum(1 for p in all_plots if p.get("topo_note"))
    print("\n%d plots from %d surveys over %d villages; geometry: %s" %
          (len(g), len(all_report), len(villages), src))
    print("topology touched %d plot(s); %d untouched plots keep the sheet's exact drawing"
          % (touched, len(g) - touched))
    print("cross-village claims: %s" % (across or "none"))
    print("written to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
