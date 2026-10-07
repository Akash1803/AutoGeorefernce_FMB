"""Place the next ring of FMB sheets onto the parcels already placed in QGIS (Manual_Georeferenced/<Village>_Village).

    python auto_georef_ring.py ring --village 35_04_221 --surveys 264,71B,71A --in-order

Per sheet: PDF cross-check (title-block area, plot numbers, neighbour numbers), rigid placement on the
placed parcels, then topology for the batch (his parcels fixed) and two independent checks (shapely on the
files, QGIS on the layers). Outputs: FMB_Vector/<village>/<survey>_parcels_auto.gpkg (red in QGIS) and
Auto_Georef/<village>/batch_NN_<date>/batch_report.csv|xlsx + run.json. Akash's layers are only read.
"""
import argparse
import datetime
import json
import sys
import time

import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

from autogeoref import paths, printed_numbers, ring, ring_io, ring_topo

TOPO_GUARD = 0.02          # topology may change a survey's area by at most 2 %


def village_name(village):
    s = paths.surveys_with_sheets(village)[0]
    props = json.loads(paths.sheet_path(village, s).read_text(encoding="utf-8"))["features"][0]["properties"]
    return props["village"]


def pdf_check(village, sheet):
    """Cross-check of the converted sheet against its PDF: title-block area, plot numbers, neighbours."""
    pdf = paths.sketch_dir(village) / ("%s.pdf" % sheet.survey)
    header = printed_numbers.read_header_area(pdf) if pdf.exists() else None
    drawn = sum(p.area for _props, p in sheet.plots)
    numbered = [str(props.get("plot_no")) for props, _p in sheet.plots if props.get("numbered")]
    review = sum(1 for props, _p in sheet.plots if props.get("review_needed"))
    unnumbered = sum(1 for props, _p in sheet.plots if not props.get("numbered") and len(sheet.plots) > 1)
    out = {"pdf_area_m2": round(header, 1) if header else None, "drawn_area_m2": round(drawn, 1),
           "area_diff_pct": round(100 * (drawn / header - 1), 1) if header else None,
           "plot_numbers": ", ".join(numbered), "plots_to_check": review, "plots_without_number": unnumbered,
           "printed_neighbours": ", ".join(sorted({L.number for L in sheet.labels if L.number != sheet.survey.upper()}))}
    notes = []
    if header and abs(drawn / header - 1) > 0.10:
        notes.append("PDF title block says %.0f m2, the drawing is %.0f m2" % (header, drawn))
    if review:
        notes.append("%d plot number(s) to confirm" % review)
    out["pdf_note"] = "; ".join(notes)
    return out


def wanted_list(a):
    return [x.strip() for x in a.surveys.split(",")] if a.surveys else []


def cmd_ring(a):
    t0 = time.time()
    vname = village_name(a.village)
    seeds, sources, seed_plots = ring_io.load_seeds(vname, include_red=a.build_on_red)
    if a.build_on_red and wanted_list(a):
        for s_ in wanted_list(a):                    # the sheets being re-placed are not their own ground
            seeds.pop(s_, None); sources.pop(s_, None)
        seed_plots = seed_plots[~seed_plots.survey_no.isin(wanted_list(a))]
    if not seeds:
        sys.exit("no placed parcels in QGIS group Manual_Georeferenced/%s" % ring_io.group_name(vname))
    print("%s: %d placed surveys in QGIS (%s)" % (vname, len(seeds), ", ".join(sorted(seeds, key=paths.survey_sort_key))))
    table = printed_numbers.village_table(a.village, str(paths.PROJECT / "Auto_Georef" / a.village / "printed_numbers.json"))
    wanted = [s.strip() for s in a.surveys.split(",")] if a.surveys else None
    if wanted and not a.no_qgis:
        print("cleared earlier red layers:", ring_io.clear_unapproved(a.village, wanted))
    decisions, sheet_objs, ctx = ring.place_ring(a.village, seeds, table, size=a.size, ring=wanted, in_order=a.in_order)
    d = ring_io.batch_dir(a.village)

    rows, placed = [], []
    for dec in decisions:
        sh = sheet_objs[dec.survey]
        row = {"survey": dec.survey, "status": dec.status, "note": dec.reason}
        row.update(dec.metrics)
        row.update(pdf_check(a.village, sh))
        if dec.status == "placed":
            row["rotation_deg"] = round(dec.theta, 2)
            placed.append(ring.plots_on_ground(sh, dec.theta, dec.t))
        rows.append(row)

    files, check1, check2 = [], None, None
    if placed:
        batch = gpd.GeoDataFrame(pd.concat(placed, ignore_index=True), geometry="geometry", crs=ring.UTM)
        batch["survey_no"] = batch["survey_no"].astype(str)
        rigid_area = batch.groupby("survey_no").geometry.apply(lambda g: unary_union(g.values).area)
        cleaned, rep = ring_topo.clean_batch(seed_plots.to_crs(ring.UTM), batch, ctx.agreed)
        cleaned_ll = ring_topo.final_pass(seed_plots, ring_topo.snap_to_seeds(cleaned.to_crs(4326), seed_plots))
        moved = rep.groupby("survey_no").moved_max_m.max()
        # guard (rule of 2026-09-21, broken by 72 on 2026-10-07: cut from 14,147 to 4,868 m2): topology may only
        # close seams; a survey it would change by more than TOPO_GUARD keeps its rigid sheet and is flagged
        rigid_ll = batch.to_crs(4326)
        carved = {}
        for sv in sorted(set(batch.survey_no)):
            after = unary_union(cleaned_ll[cleaned_ll.survey_no == sv].to_crs(ring.UTM).geometry.values).area
            if abs(after / rigid_area[sv] - 1) > TOPO_GUARD:
                carved[sv] = 100 * (after / rigid_area[sv] - 1)
        if carved:
            cleaned_ll = pd.concat([cleaned_ll[~cleaned_ll.survey_no.isin(list(carved))],
                                    rigid_ll[rigid_ll.survey_no.isin(list(carved))]], ignore_index=True)
            cleaned_ll = gpd.GeoDataFrame(cleaned_ll, geometry="geometry", crs=4326)
        for r in rows:
            if r["status"] != "placed":
                continue
            part = cleaned_ll[cleaned_ll.survey_no == r["survey"]]
            area_after = unary_union(part.to_crs(ring.UTM).geometry.values).area
            r["topology_max_move_m"] = float(moved.get(r["survey"], 0.0))
            if r["survey"] in carved:
                r["topology_max_move_m"] = 0.0
                r["note"] = ("topology would change it by %.0f %%: it overlaps its neighbours, so its position is wrong "
                             "somewhere - left as the rigid sheet for you to see; " % carved[r["survey"]]) + (r.get("note") or "")
                r["confidence"] = "conflict"
            r["area_change_by_topology_m2"] = round(area_after - rigid_area[r["survey"]], 1)
            files.append((r["survey"], ring_io.write_survey(a.village, r["survey"], part)))
        written = gpd.GeoDataFrame(pd.concat([gpd.read_file(p, layer="parcels") for _s, p in files], ignore_index=True),
                                   geometry="geometry", crs=4326)
        check1 = ring_topo.verify_files(seed_plots, written)
        print("check 1 (files, shapely):", {k: v for k, v in check1.items() if k != "overlap_pairs"})
        if not a.no_satellite:
            from shapely.geometry import box
            from autogeoref import satcheck
            polys = {s: unary_union(written[written.survey_no == s].to_crs(ring.UTM).geometry.values) for s, _p in files}
            img = satcheck.Imagery(a.village, list(polys.values()), d.name, out_dir=d)
            win = box(*img.bounds)
            base = [img.share(ring.largest(g)) for g in seeds.values() if ring.largest(g).within(win)]
            baseline = float(sorted(base)[len(base) // 2]) if base else 0.0
            for r in rows:
                if r["survey"] in polys:
                    r.update(img.verdict(ring.largest(polys[r["survey"]]), baseline))
    if files and not a.no_qgis:
        print("QGIS:", ring_io.qgis_show(vname, files))
        check2 = ring_topo.verify_qgis(vname, [s for s, _p in files])
        print("check 2 (QGIS):", {k: v for k, v in check2.items() if k != "overlap_pairs"})
    ring_io.write_report(d, rows)
    ring_io.save_run(d, {"village": a.village, "when": datetime.datetime.now(), "seeds": sources,
                         "ring": [r["survey"] for r in rows], "in_order": a.in_order,
                         "check_1_files": check1, "check_2_qgis": check2, "seconds": round(time.time() - t0)})
    print("batch folder:", d)
    for r in rows:
        print("  %-5s %-6s %-5s %s | %s | %s" % (r["survey"], r["status"], r.get("confidence", ""), r.get("support_by", ""),
                                             r.get("note") or "", r.get("pdf_note") or ""))
    print("done in %.0f s" % (time.time() - t0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("ring", help="place the next ring of sheets")
    r.add_argument("--village", required=True)
    r.add_argument("--size", type=int, default=3)
    r.add_argument("--surveys", help="place exactly these surveys (comma separated)")
    r.add_argument("--in-order", action="store_true", help="place one by one in the order given")
    r.add_argument("--build-on-red", action="store_true", help="red (not yet approved) parcels of earlier batches count as placed")
    r.add_argument("--no-qgis", action="store_true")
    r.add_argument("--no-satellite", action="store_true", help="skip the Google Satellite check")
    r.set_defaults(fn=cmd_ring)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
