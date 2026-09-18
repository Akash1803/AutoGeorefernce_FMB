"""python -m autogeoref.georef_village 35_04_077 [--raster] [--topology] [--refit 46B,47A] ..."""
import argparse
import sys

from . import engine, gcp, paths, review


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
    summary = engine.run(args.village, do_raster=args.raster, do_topology=args.topology,
                         do_review=args.review, project=args.project)
    rows = review.read_status(args.village)
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
