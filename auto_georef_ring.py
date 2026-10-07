"""Place the next ring of FMB sheets onto the parcels already placed in QGIS (Manual_Georeferenced/<Village>_Village).

    python auto_georef_ring.py ring --village 35_04_221 [--size 3] [--surveys 74,73] [--no-qgis]

Outputs: FMB_Vector/<village>/<survey>_parcels_auto.gpkg per placed sheet (shown red in QGIS for review)
and Auto_Georef/<village>/batch_NN_<date>/batch_report.csv|xlsx + run.json. Akash's layers are only read.
"""
import argparse
import datetime
import json
import sys
import time

import geopandas as gpd
import shapely
from shapely.ops import unary_union

from autogeoref import paths, printed_numbers, ring, ring_io


def village_name(village):
    s = paths.surveys_with_sheets(village)[0]
    props = json.loads(paths.sheet_path(village, s).read_text(encoding="utf-8"))["features"][0]["properties"]
    return props["village"]


def cmd_ring(a):
    t0 = time.time()
    vname = village_name(a.village)
    seeds, sources = ring_io.load_seeds(vname)
    if not seeds:
        sys.exit("no placed parcels in QGIS group Manual_Georeferenced/%s" % ring_io.group_name(vname))
    print("%s: %d placed surveys in QGIS (%s)" % (vname, len(seeds), ", ".join(sorted(seeds, key=paths.survey_sort_key))))
    table = printed_numbers.village_table(a.village, str(paths.PROJECT / "Auto_Georef" / a.village / "printed_numbers.json"))
    wanted = [s.strip() for s in a.surveys.split(",")] if a.surveys else None
    if wanted and not a.no_qgis:
        print("cleared earlier red layers:", ring_io.clear_unapproved(a.village, wanted))
    decisions, sheet_objs, ctx = ring.place_ring(a.village, seeds, table, size=a.size, ring=wanted)
    d = ring_io.batch_dir(a.village)
    seed_union = unary_union(list(seeds.values()))
    rows, files = [], []
    for dec in decisions:
        sh = sheet_objs[dec.survey]
        row = {"survey": dec.survey, "status": dec.status, "reason": dec.reason}
        row.update(dec.metrics)
        stated = sh.meta.get("stated_area_sqm")
        row["sheet_area_m2"] = round(sh.polygon.area, 1)
        row["stated_area_m2"] = round(stated, 1) if stated else None
        if dec.status == "placed":
            plots = ring.plots_on_ground(sh, dec.theta, dec.t)
            G = ring.ground_polygon(sh, dec.theta, dec.t)
            row["overlap_with_placed_m2"] = round(G.intersection(seed_union).area, 1)
            row["rotation_deg"] = round(dec.theta, 2)
            row["plots"] = len(plots)
            row["plots_to_check_number"] = int(plots.get("review_needed", gpd.GeoSeries([])).fillna(False).astype(bool).sum())
            files.append((dec.survey, ring_io.write_survey(a.village, dec.survey, plots)))
        rows.append(row)
    # options of the asked sheets, for Akash to pick from
    opt_rows = []
    for dec in decisions:
        if dec.status == "ask":
            for k, (th, t, shared, how) in enumerate(dec.options, start=1):
                opt_rows.append({"survey": dec.survey, "option": k, "label": "%s option %d" % (dec.survey, k),
                                 "shared_m": round(shared, 1), "how": how, "rotation_deg": round(th, 2),
                                 "geometry": ring.ground_polygon(sheet_objs[dec.survey], th, t)})
    if not a.no_satellite:
        from autogeoref import satcheck
        placed_polys = {r["survey"]: ring.ground_polygon(sheet_objs[r["survey"]], dec.theta, dec.t)
                        for r, dec in zip(rows, decisions) if dec.status == "placed"}
        polys = list(placed_polys.values()) + [o["geometry"] for o in opt_rows]
        if polys:
            print("satellite window ...")
            img = satcheck.Imagery(a.village, polys, d.name, out_dir=d)
            from shapely.geometry import box
            win = box(*img.bounds)
            base = [img.share(ring.largest(g)) for s, g in seeds.items() if ring.largest(g).within(win)]
            baseline = float(sorted(base)[len(base) // 2]) if base else 0.0
            for r in rows:
                if r["survey"] in placed_polys:
                    r.update(img.verdict(placed_polys[r["survey"]], baseline))
            for o in opt_rows:
                o["sat_share"] = round(img.share(o["geometry"]), 2)
                o["sat_share_his_parcels"] = round(baseline, 2)
    if opt_rows:
        og = gpd.GeoDataFrame(opt_rows, geometry="geometry", crs=32644).to_crs(4326)
        og.to_file(d / "ask_options.gpkg", layer="ask_options", driver="GPKG")
        og.drop(columns="geometry").to_csv(d / "ask_options.csv", index=False)
        if not a.no_qgis:
            print("QGIS:", ring_io.qgis_show_options(vname, d / "ask_options.gpkg", "ask_options_" + d.name.split("_2")[0]))
    ring_io.write_report(d, rows)
    ring_io.save_run(d, {"village": a.village, "when": datetime.datetime.now(), "seeds": sources,
                         "ring": [r["survey"] for r in rows], "seconds": round(time.time() - t0)})
    if files and not a.no_qgis:
        print("QGIS:", ring_io.qgis_show(vname, files))
    print("batch folder:", d)
    for r in rows:
        print("  %-5s %-6s %s" % (r["survey"], r["status"], r["reason"] or "support %s m (%s)" % (r.get("support_m"), r.get("support_by"))))
    print("done in %.0f s" % (time.time() - t0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("ring", help="place the next ring of sheets")
    r.add_argument("--village", required=True)
    r.add_argument("--size", type=int, default=3)
    r.add_argument("--surveys", help="place exactly these surveys (comma separated)")
    r.add_argument("--no-qgis", action="store_true")
    r.add_argument("--no-satellite", action="store_true", help="skip the Google Satellite check")
    r.set_defaults(fn=cmd_ring)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
