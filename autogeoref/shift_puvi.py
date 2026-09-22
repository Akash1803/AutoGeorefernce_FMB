"""Write a shift-corrected copy of the Puvi survey polygons for one or more villages.

    python -m autogeoref.shift_puvi 35_04_074 35_04_052 --buffer-only
    python -m autogeoref.shift_puvi --stretch          # Thailavaram to Thirukatchur

Control is the team's own placements: every `<survey>_parcels_modified.gpkg` in the village folder,
plus any file named with `--control`, which is how the Thailavaram merged layer comes in. The
displacement between a control parcel and its Puvi twin is measured at the centroid, the field is
fitted by `shiftfit`, and every Puvi polygon in scope is translated by it. Nothing is rotated,
scaled or reshaped, nothing existing is edited, and a village with no control is written out
unchanged and marked so.
"""
import argparse
import datetime
import logging
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import affinity
from shapely.ops import unary_union
from shapely.validation import make_valid

from . import paths, shiftfit

log = logging.getLogger(__name__)

MAX_CONTROL_M = 200.0   # a control parcel further than this from its Puvi twin is a key clash
STRETCH = ["35_04_052", "35_04_054", "35_04_056", "35_04_076", "35_04_077", "35_04_074"]
UTM = 32644


def _valid(geom):
    if geom is None or geom.is_empty:
        return None
    return geom if geom.is_valid else make_valid(geom)


def puvi_polygons(village):
    """Every Puvi survey polygon of a village, dissolved per survey, in EPSG:32644."""
    files = sorted((paths.PUVI / village[:2] / village[3:5] / village[6:] / "vector").glob("*_vector.shp"))
    if not files:
        return gpd.GeoDataFrame(columns=["key", "geometry"], geometry="geometry", crs=UTM)
    g = gpd.read_file(files[0]).to_crs(UTM)
    g["key"] = g["survey_no"].map(shiftfit.normalise)
    g["geometry"] = [_valid(x) for x in g.geometry]
    g = g[g.geometry.notna()]
    return g


def control_parcels(village, extra_files=()):
    """The team's placements for this village as {survey key: geometry}, in EPSG:32644."""
    out = {}
    for f in sorted(paths.vector_dir(village).glob("*_parcels_modified.gpkg")):
        key = shiftfit.normalise(f.name.split("_parcels_modified")[0])
        g = gpd.read_file(f)
        if g.crs is not None and g.crs.to_epsg() != UTM:
            g = g.to_crs(UTM)
        body = unary_union([x for x in (_valid(v) for v in g.geometry) if x is not None])
        if not body.is_empty:
            out[key] = body
    for path in extra_files:
        g = gpd.read_file(path)
        if g.crs is not None and g.crs.to_epsg() != UTM:
            g = g.to_crs(UTM)
        if "survey_no" not in g.columns:
            log.warning("%s has no survey_no column, skipped as control", path)
            continue
        if "village_code" in g.columns:
            # survey numbers repeat across villages: a control file only speaks for its own village
            g = g[g["village_code"].astype(str) == village]
            if not len(g):
                continue
        g["key"] = g["survey_no"].map(shiftfit.normalise)
        for key, sub in g.groupby("key"):
            if key in ("NAN", "NONE", ""):
                continue
            body = unary_union([x for x in (_valid(v) for v in sub.geometry) if x is not None])
            if not body.is_empty:
                out[key] = body
    return out


def in_buffer_keys(village):
    """Survey keys of the village that touch the 30 m rail buffer, from the team's own layer."""
    b = gpd.read_file(paths.BUFFER_GPKG, layer="vector_in_buffer_30m")
    b = b[b["village_code"].astype(str) == village]
    return {shiftfit.normalise(s) for s in b["survey_no"]}


def village_name(village, puvi):
    for col in ("village", "vill_name", "village_na"):
        if col in puvi.columns and len(puvi):
            return str(puvi[col].iloc[0])
    return village


def run_village(village, out_dir, buffer_only=True, extra_control=(), stamp=None):
    """Write the corrected polygons for one village. Returns the report row."""
    puvi = puvi_polygons(village)
    if not len(puvi):
        log.warning("%s: no Puvi vector found", village)
        return None
    control = control_parcels(village, extra_control)
    scope = in_buffer_keys(village) if buffer_only else set(puvi["key"])
    targets = puvi[puvi["key"].isin(scope)].copy()
    if not len(targets):
        log.warning("%s: no Puvi polygon is in scope", village)
        return None

    pairs, rejected = [], []
    by_key = {k: sub.geometry.iloc[0] for k, sub in puvi.groupby("key")}
    for key, hand_geom in control.items():
        if key not in by_key:
            continue
        d = float(np.linalg.norm(np.array(hand_geom.centroid.coords[0])
                                 - np.array(by_key[key].centroid.coords[0])))
        if d > MAX_CONTROL_M:
            # the same survey number in two villages, or a mis-placed parcel: not control
            rejected.append((key, round(d)))
            continue
        pairs.append((key, hand_geom))
    if rejected:
        log.warning("%s: %d control parcel(s) ignored, more than %.0f m from their Puvi twin: %s",
                    village, len(rejected), MAX_CONTROL_M, sorted(rejected, key=lambda r: -r[1])[:5])
    cpoints, cdisp = [], []
    for key, hand_geom in pairs:
        p = by_key[key]
        cpoints.append(np.array(p.centroid.coords[0]))
        cdisp.append(np.array(hand_geom.centroid.coords[0]) - np.array(p.centroid.coords[0]))
    cpoints = np.array(cpoints) if cpoints else np.zeros((0, 2))
    cdisp = np.array(cdisp) if cdisp else np.zeros((0, 2))

    tpoints = np.array([g.centroid.coords[0] for g in targets.geometry])
    shifts, method = shiftfit.fit(cpoints, cdisp, tpoints)
    acc = shiftfit.accuracy(cpoints, cdisp)

    moved, source = [], []
    control_keys = {k for k, _ in pairs}
    for (i, row), s in zip(targets.iterrows(), shifts):
        if row["key"] in control_keys:
            moved.append(control[row["key"]])          # the team placed it: their geometry stands
            source.append("team placement")
        else:
            moved.append(affinity.translate(row.geometry, float(s[0]), float(s[1])))
            source.append("puvi shifted")
    out = targets.copy()
    out["geometry"] = moved
    out["village_code"] = village
    out["survey_no"] = out["key"]
    out["shift_x_m"] = [round(float(s[0]), 2) for s in shifts]
    out["shift_y_m"] = [round(float(s[1]), 2) for s in shifts]
    out["shift_m"] = [round(float(np.hypot(*s)), 2) for s in shifts]
    out["source"] = source
    out["fit_method"] = method
    out["control_parcels"] = len(cpoints)
    out["expected_error_m"] = acc["after_median_m"] if acc["after_median_m"] is not None else ""
    out["puvi_error_before_m"] = acc["before_median_m"] if acc["before_median_m"] is not None else ""
    out["run_at"] = stamp or shiftfit.run_stamp()

    keep = ["village_code", "survey_no", "source", "fit_method", "control_parcels",
            "shift_x_m", "shift_y_m", "shift_m", "expected_error_m", "puvi_error_before_m",
            "run_at", "geometry"]
    out = out[[c for c in keep if c in out.columns]]
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / ("%s_puvi_shifted.geojson" % village)
    out.to_crs(4326).to_file(dest, driver="GeoJSON", COORDINATE_PRECISION=8)
    log.info("%s: %d polygons -> %s (%s, %d control)", village, len(out), dest, method, len(cpoints))

    return {"village_code": village, "village_name": village_name(village, puvi),
            "control": len(cpoints), "targets": len(out), "method": method,
            "mean_shift_x_m": round(float(cdisp.mean(0)[0]), 2) if len(cdisp) else "",
            "mean_shift_y_m": round(float(cdisp.mean(0)[1]), 2) if len(cdisp) else "",
            "run_at": out["run_at"].iloc[0], **{k: v for k, v in acc.items() if k != "control"}}


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="shift_puvi", description=__doc__.splitlines()[0])
    ap.add_argument("villages", nargs="*", help="village codes, e.g. 35_04_074")
    ap.add_argument("--stretch", action="store_true", help="the Thailavaram to Thirukatchur villages")
    ap.add_argument("--buffer-only", action="store_true", default=True,
                    help="only surveys touching the 30 m rail buffer (default)")
    ap.add_argument("--whole-village", dest="buffer_only", action="store_false",
                    help="every survey of the village, not just the buffer")
    ap.add_argument("--control", action="append", default=[],
                    help="extra control file with a survey_no column (repeatable)")
    ap.add_argument("--out", default="", help="output folder (default: <project>\\Puvi_Shifted\\<date>)")
    args = ap.parse_args(argv)

    villages = list(args.villages) + (STRETCH if args.stretch else [])
    if not villages:
        ap.error("name at least one village, or pass --stretch")
    stamp = shiftfit.run_stamp()
    out_dir = __import__("pathlib").Path(args.out) if args.out else shiftfit.output_dir()
    rows = []
    for v in dict.fromkeys(villages):
        row = run_village(v, out_dir, buffer_only=args.buffer_only,
                          extra_control=args.control, stamp=stamp)
        if row:
            rows.append(row)
    if not rows:
        log.error("nothing written")
        return 1
    merged = pd.concat([gpd.read_file(out_dir / ("%s_puvi_shifted.geojson" % r["village_code"]))
                        for r in rows], ignore_index=True)
    gpd.GeoDataFrame(merged, crs=4326).to_file(out_dir / "puvi_shifted_all.geojson",
                                               driver="GeoJSON", COORDINATE_PRECISION=8)
    shiftfit.write_report(rows, out_dir / "shift_report.csv")
    print("\n%-16s %8s %8s %-12s %10s %10s" % ("village", "control", "parcels", "method", "before", "after"))
    for r in rows:
        print("%-16s %8d %8d %-12s %9s m %9s m"
              % (r["village_name"][:16], r["control"], r["targets"], r["method"],
                 r["before_median_m"] if r["before_median_m"] is not None else "-",
                 r["after_median_m"] if r["after_median_m"] is not None else "-"))
    print("\nwritten to %s" % out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
