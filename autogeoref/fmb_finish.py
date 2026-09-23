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
CONFORM_STEP_M = 3.0    # boundary correspondence sample spacing
CONFORM_K = 8           # neighbours in the conform field
CONFORM_SMOOTH_M = 5.0  # keeps the field smooth between correspondence points
AREA_HARD = (1.0 / 3.0, 3.0)   # outside the 3x gate the sheet and the base are different ground


def conform_survey(rows, base_geom):
    """Warp one placed sheet so its outline IS its base parcel's outline, plots carried inside.

    Akash, 2026-09-23: both layers should fit. The base outline is authoritative; the sheet
    contributes the interior subdivision. Every distinct vertex moves once (the node table), so
    shared plot lines stay one line; the exact fit at the end clips to the base parcel and gives
    whatever remains of it to the plot bordering it most, so the union equals the base parcel.
    Printed sheet lengths no longer survive exactly; the rigid pose is what preserved them, and
    shape_match in the report still says how much the two drawings disagree.
    Returns how far the outline had to move, in metres.
    """
    geoms = [fmb_on_base._valid(r["geometry"]) for r in rows]
    body = fmb_on_base._valid(unary_union([g for g in geoms if g is not None]))
    if body is None or body.is_empty:
        return 0.0
    bb = base_geom.boundary
    pts, disp = [], []
    for poly in shift_puvi._polygons_of(body):
        ring = poly.exterior
        n = max(8, int(ring.length / CONFORM_STEP_M))
        for i in range(n):
            s = ring.interpolate(i * ring.length / n)
            q = bb.interpolate(bb.project(s))
            pts.append((s.x, s.y))
            disp.append((q.x - s.x, q.y - s.y))
    pts, disp = np.asarray(pts), np.asarray(disp)
    moved = float(np.linalg.norm(disp, axis=1).max()) if len(pts) else 0.0
    table = shift_puvi._warp_table(shift_puvi._nodes([g for g in geoms if g is not None]),
                                   pts, disp, k=CONFORM_K, smooth=CONFORM_SMOOTH_M)
    for r, g in zip(rows, geoms):
        if g is not None:
            r["geometry"] = fmb_on_base._valid(shift_puvi._warp_with(g, table))
    parts = [(r, r["geometry"]) for r in rows if r["geometry"] is not None]
    new_body = fmb_on_base._valid(unary_union([g for _r, g in parts]))
    if parts and new_body is not None and not new_body.is_empty:
        for (r, _g), (_r2, ng, note) in zip(parts,
                                            topology.apply_to_parts(parts, base_geom, new_body)):
            r["geometry"] = ng
            if note.startswith("kept"):
                r["topo_note"] = note
    return moved


def _adopt_sheet(village, survey, base_geom, proto, rep):
    """Bring a disagreeing sheet's plots in anyway, for conforming (Akash, 2026-09-23).

    The shape gate stays only as the 3x area rule: a sheet describing wholly different ground
    (569B's 0.98 acres against a 259 acre tank) is not stretched over it. Everything nearer than
    that is placed and conformed; the flag keeps saying how much the drawings disagreed.
    """
    sheet_gdf, outline = fmb_on_base.sheet_of(village, survey)
    if outline is None:
        return None
    pose = fmb_on_base.place(outline, base_geom)
    if pose is None or not (AREA_HARD[0] <= pose["area_ratio"] <= AREA_HARD[1]):
        return None
    new = []
    for _i, srow in sheet_gdf.iterrows():
        g = fmb_on_base.apply_pose(fmb_on_base._valid(srow.geometry), pose, base_geom)
        if g is None or g.is_empty:
            continue
        r = {k: proto[k] for k in ("village_code", "village", "survey_no", "evidence", "run_at")
             if k in proto}
        r.update({"geometry_source": "fmb sheet", "plot_no": str(srow.get("plot_no", "")),
                  "shape_match": round(pose["match"], 3),
                  "area_ratio": round(pose["area_ratio"], 3), "gcp_points": 0, "geometry": g})
        new.append(r)
    if not new:
        return None
    if rep is not None:
        rep["geometry_source"] = "fmb sheet"
        rep["plots_expected"] = len(sheet_gdf)
        rep["plots_kept"] = len(new)
        rep["flag"] = (str(rep.get("flag", "")) + "; conformed despite the disagreement").strip("; ")
    return new


def _unit_names(village, survey, base_keys):
    """Letter units of a base survey on disk: 2A and 2B for a base parcel the Puvi layer calls 2.

    The portal subdivides further than Puvi. A unit that is itself a parcel of the base
    (52A in Thirukatchur) is never treated as a piece of another survey.
    """
    import re
    out = []
    vd = paths.vector_dir(village)
    if not vd.exists():
        return out
    for f in sorted(vd.glob("%s*_parcels.geojson" % survey)):
        stem = f.name[:-len("_parcels.geojson")]
        m = re.match(r"^%s([A-Z]+)$" % re.escape(survey), stem)
        if m and stem not in base_keys:
            out.append(stem)
    return out


def _slabs(base_geom, fractions):
    """Cut the base parcel into slabs across its principal axis, one per unit, sized by area."""
    from shapely.geometry import Polygon as _P
    hull = unary_union(shift_puvi._polygons_of(base_geom)).convex_hull
    coords = np.array(hull.exterior.coords)[:, :2]
    c = coords.mean(0)
    _u, _s, vt = np.linalg.svd(coords - c)
    axis, perp = vt[0], np.array([-vt[0][1], vt[0][0]])
    tvals = (coords - c) @ axis
    t0, t1 = float(tvals.min()), float(tvals.max())
    cuts = [t0 + f * (t1 - t0) for f in fractions]
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        rect = _P([c + a * axis + 1e5 * perp, c + b * axis + 1e5 * perp,
                   c + b * axis - 1e5 * perp, c + a * axis - 1e5 * perp])
        out.append(fmb_on_base._valid(base_geom.intersection(rect)))
    return out


def _clip_within(rows):
    """No plot of an assembled survey claims another's ground; the larger one yields."""
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            a, b = rows[i]["geometry"], rows[j]["geometry"]
            if a is None or b is None:
                continue
            inter = a.intersection(b)
            if inter.is_empty or inter.area <= 0.05:
                continue
            big, small = (i, j) if a.area >= b.area else (j, i)
            if rows[big]["geometry"].area > 0 and inter.area / rows[big]["geometry"].area > topology.GUARD_FRACTION:
                continue
            cut = fmb_on_base._valid(rows[big]["geometry"].difference(rows[small]["geometry"]))
            if cut is not None and not cut.is_empty:
                rows[big]["geometry"] = cut


def _adopt_units(village, survey, base_geom, proto, rep, base_keys):
    """Assemble a base parcel from its portal units (2A and 2B making Puvi's 2).

    The units partition the parcel, so each is fitted into its own slab of the parcel across
    the principal axis, sized by drawn area, in letter order - and in reverse letter order,
    whichever the drawings overlap better. The conform pass then makes them tile it exactly.
    """
    units = []
    for u in _unit_names(village, survey, base_keys):
        gdf, outline = fmb_on_base.sheet_of(village, u)
        if outline is not None:
            units.append((u, gdf, outline))
    if not units:
        return None
    total = sum(o.area for _u, _g, o in units)
    ratio = total / base_geom.area if base_geom.area else 0.0
    if not (AREA_HARD[0] <= ratio <= AREA_HARD[1]):
        return None

    best = None
    for order in (units, units[::-1]):
        fractions = [0.0]
        for _u, _g, o in order:
            fractions.append(fractions[-1] + o.area / total)
        slabs = _slabs(base_geom, fractions)
        placed, score = [], 0.0
        for (u, gdf, outline), slab in zip(order, slabs):
            if slab is None or slab.is_empty:
                continue
            pose = fmb_on_base.place(outline, slab)
            if pose is None:
                continue
            placed.append((u, gdf, outline, slab, pose))
            score += pose["match"]
        if placed and (best is None or score > best[1]):
            best = (placed, score)
    if best is None:
        return None

    new, used = [], []
    for u, gdf, _outline, slab, pose in best[0]:
        kept = 0
        for _i, srow in gdf.iterrows():
            g = fmb_on_base.apply_pose(fmb_on_base._valid(srow.geometry), pose, slab)
            if g is None or g.is_empty:
                continue
            r = {k: proto[k] for k in ("village_code", "village", "survey_no", "evidence",
                                       "run_at") if k in proto}
            r.update({"geometry_source": "fmb sheet",
                      "plot_no": "%s/%s" % (u, srow.get("plot_no", "")),
                      "shape_match": round(pose["match"], 3), "area_ratio": round(ratio, 3),
                      "gcp_points": 0, "geometry": g})
            new.append(r)
            kept += 1
        if kept:
            used.append(u)
    if not new:
        return None
    _clip_within(new)
    if rep is not None:
        rep["geometry_source"] = "fmb sheet"
        rep["plots_expected"] = sum(len(g) for _u, g, _o in units)
        rep["plots_kept"] = len(new)
        rep["flag"] = (str(rep.get("flag", "") or "")
                       + "; assembled from units %s" % "+".join(used)).strip("; ")
    return new


def conform_village(village, plots, report, base_sub):
    """Conform every survey with a sheet to its base parcel. Returns (plots, {survey: m})."""
    base_of = {}
    for _, row in base_sub.iterrows():
        g = fmb_on_base._valid(row.geometry)
        if g is not None:
            base_of.setdefault(str(row["survey_no"]), []).append(g)
    by_survey = {}
    for p in plots:
        by_survey.setdefault(str(p["survey_no"]), []).append(p)
    rep_of = {str(r["survey_no"]): r for r in report}
    base_keys = set(base_of)
    out_plots, stats = [], {}
    for s, rows in sorted(by_survey.items()):
        src = rows[0].get("geometry_source")
        blist = base_of.get(s)
        base_geom = fmb_on_base._valid(unary_union(blist)) if blist else None
        if base_geom is None or base_geom.is_empty:
            out_plots += rows
            continue
        if src == "your placement":
            # Akash, 2026-09-23: ALL parcels sit exactly on the base, his 80 included. His own
            # files and layers stay untouched; only this deliverable refits. The sheet's plots
            # come in where a sheet exists, else his body simply takes the base parcel's shape.
            adopted = (_adopt_sheet(village, s, base_geom, rows[0], rep_of.get(s))
                       or _adopt_units(village, s, base_geom, rows[0], rep_of.get(s), base_keys))
            if adopted is not None:
                rows = adopted
            else:
                for r in rows:
                    r["geometry_source"] = "puvi base"
            rep = rep_of.get(s)
            if rep is not None:
                if adopted is None:
                    rep["geometry_source"] = "puvi base"
                rep["flag"] = (str(rep.get("flag", "") or "")
                               + "; was your placement, refitted to the base").strip("; ")
        elif src == "puvi base":
            adopted = (_adopt_sheet(village, s, base_geom, rows[0], rep_of.get(s))
                       or _adopt_units(village, s, base_geom, rows[0], rep_of.get(s), base_keys))
            if adopted is None:
                out_plots += rows          # no sheet, or truly different ground: the base stands
                continue
            rows = adopted
        stats[s] = round(conform_survey(rows, base_geom), 2)
        out_plots += rows
    for r in report:
        r["conform_max_m"] = stats.get(str(r["survey_no"]), "")
    for p in out_plots:
        p["conform_max_m"] = stats.get(str(p["survey_no"]), "")
    return out_plots, stats


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


def settle_cross_village(all_plots, topo_rows):
    """Clip the land two villages both claim, plot by plot; a village with Akash's own
    placements keeps its ground, otherwise the larger body yields. Guarded like every clip."""
    by_village, hand_villages = {}, set()
    for pl in all_plots:
        v = str(pl["village_code"])
        by_village.setdefault(v, []).append(pl)
        if pl.get("geometry_source") == "your placement":
            hand_villages.add(v)
    bodies = {v: fmb_on_base._valid(unary_union([pl["geometry"] for pl in rows]))
              for v, rows in by_village.items()}
    for a in sorted(by_village):
        for b in sorted(by_village):
            if b <= a or bodies[a] is None or bodies[b] is None:
                continue
            inter = bodies[a].intersection(bodies[b])
            if inter.is_empty or inter.area <= 1.0:
                continue
            if a in hand_villages and b not in hand_villages:
                yielder, keeper = b, a
            elif b in hand_villages and a not in hand_villages:
                yielder, keeper = a, b
            else:
                yielder, keeper = (a, b) if bodies[a].area >= bodies[b].area else (b, a)
            keep_body = bodies[keeper]
            by_survey = {}
            for pl in by_village[yielder]:
                by_survey.setdefault(str(pl["survey_no"]), []).append(pl)
            for s, rows in sorted(by_survey.items()):
                if rows[0].get("geometry_source") == "your placement":
                    continue
                body_s = fmb_on_base._valid(unary_union([r["geometry"] for r in rows]))
                if body_s is None:
                    continue
                piece = body_s.intersection(keep_body)
                if piece.is_empty or piece.area <= 1.0:
                    continue
                if body_s.area > 0 and piece.area / body_s.area > topology.GUARD_FRACTION:
                    topo_rows.append({"village_code": yielder, "survey_no": s,
                                      "action": "refused",
                                      "against": "cross-village clip would take %.0f %%"
                                                 % (100 * piece.area / body_s.area),
                                      "sqm": round(piece.area, 1)})
                    continue
                new_body = fmb_on_base._valid(body_s.difference(keep_body))
                if new_body is None or new_body.is_empty:
                    continue
                parts = [(r, r["geometry"]) for r in rows]
                for (r, _g), (_r2, ng, note) in zip(parts,
                                                    topology.apply_to_parts(parts, new_body, body_s)):
                    r["geometry"] = ng
                    if note:
                        r["topo_note"] = (str(r.get("topo_note", "")) + "; " + note).strip("; ")
                topo_rows.append({"village_code": yielder, "survey_no": s,
                                  "action": "cross-village clip", "against": keeper,
                                  "sqm": round(piece.area, 1)})
            bodies[yielder] = fmb_on_base._valid(
                unary_union([pl["geometry"] for pl in by_village[yielder]]))
    return all_plots


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
    """Valid in the written CRS, not only in metres: snap, then repair what remains.

    The snap grid matches COORDINATE_PRECISION=8: snapping finer than the writer rounds
    left 4 dissolved surveys invalid on 2026-09-23.
    """
    g = shift_puvi._snap_or_keep(geom, grid=1e-8)
    if g is not None and not g.is_empty and not g.is_valid:
        g = g.buffer(0)
    return g


REPORT_COLUMNS = fmb_on_base.REPORT_COLUMNS + ["topo_note", "conform_max_m"]


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
        plots, stats = conform_village(v, plots, report, sub)
        if stats:
            log.info("%s: %d survey(s) conformed to the base outline, median %.2f m, max %.2f m",
                     v, len(stats), float(np.median(list(stats.values()))), max(stats.values()))
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
    settle_cross_village(all_plots, topo_rows)
    across = cross_village_overlap(all_plots)
    for (a, b), area in sorted(across.items()):
        log.info("village overlap %s / %s: %.1f m2 still standing after the settle", a, b, area)

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
    placed_now = {(str(r["village_code"]), str(r["survey_no"])) for r in all_report
                  if r.get("geometry_source") == "fmb sheet"}
    all_unplaced = [u for u in all_unplaced
                    if (str(u["village_code"]), str(u["survey_no"])) not in placed_now]
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
