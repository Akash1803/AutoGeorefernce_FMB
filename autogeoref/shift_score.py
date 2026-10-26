"""Score a corrected Puvi layer against the team's own placements.

    python -m autogeoref.shift_score <layer.geojson> [<another.geojson> ...]

Every figure is measured on the parcels where a hand placement exists, which is the only truth
there is. The accuracy floor is the team's own placements disagreeing with each other by 1 to 3 m,
so nothing here can be read finer than that.
"""
import argparse
import logging
import sys

import geopandas as gpd
import numpy as np
from shapely.validation import make_valid

from . import shift_puvi

log = logging.getLogger(__name__)

UTM = 32644
THAILAVARAM_MERGED = r"D:\Projects\Tambaram_Chengalpattu\Thailavaram_merged fmb.geojson"
SCORED_VILLAGES = ("35_04_052", "35_04_077", "35_04_074")


def _valid(geom):
    if geom is None or geom.is_empty:
        return None
    return geom if geom.is_valid else make_valid(geom)


def rmse(values):
    values = np.asarray(values, float)
    return float(np.sqrt(np.mean(np.square(values)))) if len(values) else float("nan")


def truth_for(village, control_files=()):
    """{survey: (hand centroid, puvi centroid)} for every parcel the team has placed."""
    hand = {k: _valid(v) for k, v in shift_puvi.control_parcels(village, control_files).items()}
    puvi = {k: _valid(s.geometry.iloc[0])
            for k, s in shift_puvi.puvi_polygons(village).groupby("key")}
    out = {}
    for k, h in hand.items():
        p = puvi.get(k)
        if h is None or p is None:
            continue
        t = np.array(h.centroid.coords[0])
        s = np.array(p.centroid.coords[0])
        if np.linalg.norm(t - s) > shift_puvi.MAX_CONTROL_M:
            continue          # the same survey number in two villages, not a placement error
        out[k] = (t, s)
    return out


def score_layer(path, villages=SCORED_VILLAGES, control_files=(THAILAVARAM_MERGED,)):
    """Error of each scored parcel in `path`, and of raw Puvi, village by village."""
    g = gpd.read_file(path).to_crs(UTM)
    rows = []
    for v in villages:
        truth = truth_for(v, control_files if v == "35_04_052" else ())
        if not truth:
            continue
        placed = {}
        for _, r in g[g["village_code"].astype(str) == v].iterrows():
            geom = _valid(r.geometry)
            if geom is not None:
                placed[str(r["survey_no"])] = np.array(geom.centroid.coords[0])
        for k, (t, s) in truth.items():
            if k not in placed:
                continue
            rows.append({"village": v, "survey": k,
                         "err_puvi_m": float(np.linalg.norm(t - s)),
                         "err_layer_m": float(np.linalg.norm(t - placed[k]))})
    return rows


def summarise(rows):
    puvi = np.array([r["err_puvi_m"] for r in rows])
    layer = np.array([r["err_layer_m"] for r in rows])
    return {"parcels": len(rows), "puvi_rmse_m": round(rmse(puvi), 2),
            "layer_rmse_m": round(rmse(layer), 2),
            "puvi_median_m": round(float(np.median(puvi)), 2) if len(puvi) else None,
            "layer_median_m": round(float(np.median(layer)), 2) if len(layer) else None,
            "better": int((layer < puvi - 0.2).sum()), "worse": int((layer > puvi + 0.2).sum()),
            "within_3m": int((layer <= 3).sum()), "within_10m": int((layer <= 10).sum())}


def rail_summary(path):
    """How the railway land sits across the track, per village, in this layer."""
    g = gpd.read_file(path).to_crs(UTM)
    line = shift_puvi._rail_line()
    out = {}
    for v in sorted(set(g["village_code"].astype(str))):
        strips = shift_puvi.rail_strips(g[g["village_code"].astype(str) == v], v, line)
        if strips:
            out[v] = round(float(np.median([s["centre"] for s in strips])), 1)
    return out


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="shift_score", description=__doc__.splitlines()[0])
    ap.add_argument("layers", nargs="+", help="corrected parcel layers to score")
    ap.add_argument("--rail", action="store_true", help="also report the railway land offsets")
    args = ap.parse_args(argv)

    print("\n%-44s %8s %10s %10s %9s %8s" % ("layer", "parcels", "RMSE", "median", "better", "worse"))
    for path in args.layers:
        rows = score_layer(path)
        if not rows:
            print("%-44s  nothing could be scored" % path[-44:])
            continue
        s = summarise(rows)
        print("%-44s %8d %8.2f m %8.2f m %9d %8d"
              % (path[-44:], s["parcels"], s["layer_rmse_m"], s["layer_median_m"], s["better"], s["worse"]))
        for v in SCORED_VILLAGES:
            sub = [r for r in rows if r["village"] == v]
            if sub:
                t = summarise(sub)
                print("    %-40s %8d %8.2f m %8.2f m %9d %8d"
                      % (v, t["parcels"], t["layer_rmse_m"], t["layer_median_m"], t["better"], t["worse"]))
        if args.rail:
            print("    railway land across the track: %s" % rail_summary(path))
    first = score_layer(args.layers[0])
    if first:
        s = summarise(first)
        print("\nPuvi as it is, on the same %d parcels: RMSE %.2f m, median %.2f m"
              % (s["parcels"], s["puvi_rmse_m"], s["puvi_median_m"]))
        print("The floor is 1 to 3 m: that is how much the team's own placements disagree with each other.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
