"""python -m autogeoref.georef_village 35_04_077 [--raster] [--topology] [--refit 46B,47A] ..."""
import argparse
import sys

import datetime
import logging

from . import config as configmod, engine, evalrows, gcp, paths, raster, references as refmod, report, review, transcribe

log = logging.getLogger(__name__)


def parse(argv=None):
    ap = argparse.ArgumentParser(prog="georef_village", description="Auto-georeference one village.")
    ap.add_argument("village", help="village code, e.g. 35_04_077")
    ap.add_argument("--raster", action="store_true", help="re-export and pin the satellite GeoTIFF")
    ap.add_argument("--topology", action="store_true",
                    help="run the topology fix on the placed parcels")
    ap.add_argument("--refit", metavar="SURVEYS",
                    help="comma list: re-fit from the team's edited .points")
    ap.add_argument("--reset-gcp", metavar="SURVEY", dest="reset_gcp",
                    help="drop the team GCP rows for one survey")
    ap.add_argument("--review", dest="review", action="store_true", default=True)
    ap.add_argument("--no-review", dest="review", action="store_false",
                    help="skip the QGIS review group")
    ap.add_argument("--project", help="path to the .qgz the review group belongs in")
    ap.add_argument("--tracker", action="store_true",
                    help="also add the Auto colour/note/run columns to the tracker workbook")
    ap.add_argument("--config", help="JSON config (default configs/default.json or $AUTOGEOREF_CONFIG)")
    ap.add_argument("--rows", dest="rows", action="store_true", default=True,
                    help="append one feature row per placed parcel to _logs/eval/rows.csv (default)")
    ap.add_argument("--no-rows", dest="rows", action="store_false", help="skip the evaluation rows")
    ap.add_argument("--render-sheets", dest="render_sheets", action="store_true",
                    help="render every sheet's drawing area to the sheet-render cache for the readers, and stop")
    ap.add_argument("--transcribe", action="store_true",
                    help="merge the two readers' readings into the neighbour table, review CSV, agreement report and seed plan, and stop")
    ap.add_argument("--report", action="store_true",
                    help="write the evaluation report from every row logged so far and stop")
    ap.add_argument("--purge-cache", dest="purge_cache", action="store_true",
                    help="delete the imagery cache (windows and log) and stop")
    ap.add_argument("--migrate-cache", dest="migrate_cache", action="store_true",
                    help="move a pre-PR0 satellite window into the cache and repoint its pin")
    return ap.parse_args(argv)


def _corridor_rows(current_village, current_rows):
    """Every village's latest verdicts, so the corridor worklist stays whole."""
    out = {current_village: current_rows}
    root = paths.PROJECT / "FMB_Vector"
    if root.exists():
        for d in sorted(root.iterdir()):
            if d.is_dir() and d.name != current_village and (d / "georef_status.csv").exists():
                out[d.name] = review.read_status(d.name)
    return out


def main(argv=None):
    args = parse(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    cfg = configmod.load(args.config)
    if args.purge_cache:
        print("imagery cache removed:", raster.purge_cache(cfg))
        return 0
    if args.render_sheets:
        done = transcribe.render_village(args.village, cfg)
        print("rendered %d sheet crops under %s" % (len(done), transcribe.render_dir(cfg) / args.village))
        return 0
    if args.transcribe:
        out = transcribe.run(args.village, cfg)
        t = out["totals"]
        print("%s: %d sheets, number agreement (Jaccard) %s, side agreement %s, %d disagreements -> %s"
              % (args.village, out["sheets"], t["number_jaccard"], t["side_agreement"], out["disagreements"], out["review"]))
        print("  table:", out["override"]); print("  report:", out["report"]); print("  seeds:", out["seeds"])
        fs = out["plan"]["from_scratch"]
        print("  seeds from scratch: %d (%d already placed, to place: %s); labelled parcels: %d"
              % (fs["n"], len(fs["already_placed"]), ", ".join(fs["to_place"]) or "-", len(fs["labelled_parcels"])))
        if out["missing"]:
            print("  sheets missing from a reading:", ", ".join(out["missing"]))
        return 0
    if args.report:
        rows = evalrows.load_rows()
        out = paths.logs_dir() / "eval" / ("report_%s.md" % datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
        print("report:", report.write_report(rows, out, "all logged rows", cfg.acceptance_m))
        return 0
    if not paths.vector_dir(args.village).exists() or not paths.surveys_with_sheets(args.village):
        print("no sheets found for village %s under %s"
              % (args.village, paths.vector_dir(args.village)))
        return 2
    if args.reset_gcp:
        p = paths.gcp_csv(args.village, args.reset_gcp)
        if p.exists():
            p.unlink()
        print("team GCP record cleared for", args.reset_gcp)
    if args.refit:
        for survey in [s.strip() for s in args.refit.split(",") if s.strip()]:
            got = gcp.refit(args.village, survey)
            if got is None:
                print("%-6s no team GCPs newer than the tool's record (or fewer than two)" % survey)
                continue
            theta, t, rms, n = got
            row = {"survey": survey, "method": "team-gcp", "colour": "green", "confidence": 95,
                   "notes": "refit from %d team GCPs, rms %.2f m" % (n, rms)}
            engine.write_parcels(args.village, survey, theta, t, row)
            print("%-6s refitted from %d team GCPs: heading %.2f, rms %.2f m"
                  % (survey, n, theta, rms))
        return 0
    if args.migrate_cache:
        print("migrated:", raster.migrate_to_cache(args.village, cfg))
    summary = engine.run(args.village, do_raster=args.raster, do_topology=args.topology,
                         do_review=args.review, project=args.project)
    rows = review.read_status(args.village)
    if args.rows:
        refs = refmod.references(args.village, cfg)
        refmod.write_exclusions(args.village, refs, cfg)
        pin = raster.pinned(args.village) or {}
        feature_rows = evalrows.rows_from_run(args.village, "full_%s_%s" % (args.village, summary["run"]), "full", [],
                                              cfg, refs, imagery={"source": pin.get("source", ""),
                                                                  "sha256": pin.get("sha256", "")})
        print("  rows:", evalrows.append_rows(feature_rows), "(%d rows)" % len(feature_rows))
    print("%s: %d anchors frozen, %d placed, %d waiting, %d anchor conflicts"
          % (summary["village"], summary["anchors"], summary["placed"], summary["waiting"],
             summary["conflicts"]))
    for colour in ("green", "amber", "red"):
        names = [r["survey"] for r in rows if r.get("colour") == colour]
        if names:
            print("  %-5s %s" % (colour, ", ".join(names)))
    if args.tracker:
        print("  tracker:", review.update_tracker(
            [dict(r, village_code=args.village, run=r.get("run", "")) for r in rows]))
    print("  worklist:", review.write_worklist(_corridor_rows(args.village, rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
