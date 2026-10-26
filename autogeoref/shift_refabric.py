"""Rebuild the fine-tuned Puvi layer as ONE warp of the raw fabric, so no gaps can form.

    python -m autogeoref.shift_refabric --out "<new file>"

Akash, 2026-09-23: raw Puvi has no parcel-to-parcel gaps inside a village, but the fine-tuned
layer does. They came from the stages that moved parcels individually (the railway land placed
bodily on the track); a parcel that moves alone parts from its neighbours.

The cure keeps both properties. The fine-tuned layer says where each parcel should BE; the raw
buffer fabric says how parcels JOIN. So one displacement field per village is learned from every
parcel's fine-tuned move (its centroid displacement, no robust down-weighting: the rail strips'
different opinion is the point, not an outlier), and the RAW fabric's nodes are warped by it,
each distinct vertex once. Shared boundaries stay one line; gaps are impossible by construction.
The price is that a bodily-moved parcel lands only approximately where the fine-tune put it,
graded smoothly into its neighbours; the scorer says what that costs.
"""
import argparse
import logging
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from . import paths, shift_puvi

log = logging.getLogger(__name__)

UTM = 32644
K = 8
SMOOTH_M = 15.0
MAX_PAIR_M = 200.0     # a raw and fine row further apart than this are not the same parcel


def pair_rows(raw_v, fine_v):
    """Match each raw parcel to its fine-tuned twin: same survey number, nearest centroid.

    Never by survey number alone: Thailavaram carries 11 duplicated keys in the buffer, and
    matching them blindly invented 347 m moves once.
    """
    out = []
    fine_groups = {k: list(v) for k, v in fine_v.groupby(fine_v["survey_no"].astype(str)).groups.items()}
    for key, ridx in raw_v.groupby(raw_v["survey_no"].astype(str)).groups.items():
        fidx = fine_groups.get(key, [])
        if not fidx:
            continue
        rc = {i: np.array(raw_v.geometry[i].centroid.coords[0]) for i in list(ridx)}
        fc = {j: np.array(fine_v.geometry[j].centroid.coords[0]) for j in fidx}
        ranked = sorted((float(np.linalg.norm(rc[i] - fc[j])), i, j) for i in rc for j in fc)
        used_r, used_f = set(), set()
        for d, i, j in ranked:
            if i in used_r or j in used_f or d > MAX_PAIR_M:
                continue
            used_r.add(i)
            used_f.add(j)
            out.append((i, j))
    return out


def refabric_village(raw_v, fine_v, k=K, smooth=SMOOTH_M):
    """Warp the raw fabric of one village onto the fine-tuned positions. Returns a GeoDataFrame."""
    pairs = pair_rows(raw_v, fine_v)
    if not pairs:
        return fine_v.copy()
    pts = np.array([raw_v.geometry[i].centroid.coords[0] for i, _j in pairs])
    disp = np.array([np.array(fine_v.geometry[j].centroid.coords[0])
                     - np.array(raw_v.geometry[i].centroid.coords[0]) for i, j in pairs])
    geoms = [shift_puvi._valid(g) for g in raw_v.geometry]
    table = shift_puvi._warp_table(shift_puvi._nodes([g for g in geoms if g is not None]),
                                   pts, disp, k=k, smooth=smooth, rounds=1)
    rows = []
    matched_f = {j for _i, j in pairs}
    raw_geom = {i: g for i, g in zip(raw_v.index, geoms)}
    for i, j in pairs:
        g = raw_geom.get(i)
        row = fine_v.loc[j].to_dict()
        row["geometry"] = shift_puvi._warp_with(g, table) if g is not None else fine_v.geometry[j]
        rows.append(row)
    for j in fine_v.index:            # a fine row with no raw twin keeps its fine geometry
        if j not in matched_f:
            rows.append(fine_v.loc[j].to_dict())
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=fine_v.crs)


def close_raw_neighbour_gaps(gdf, raw_v, pairs, max_gap=8.0, rounds=2):
    """Fill any gap between parcels that TOUCH in the raw fabric; the sliver goes to the
    side that borders it most. Raw non-neighbours (real lanes and streets) are never touched."""
    from shapely.ops import unary_union as _uu
    from . import topology as _topo
    raw_geom = {i: shift_puvi._valid(g) for i, g in zip(raw_v.index, raw_v.geometry)}
    pos = {i: k for k, (i, _j) in enumerate(pairs)}
    geoms = list(gdf.geometry)                  # mutate a list; .iloc writes on a slice are lost
    raw_touch = {}
    idx = list(pos)
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if raw_geom[i] is not None and raw_geom[j] is not None                     and raw_geom[i].distance(raw_geom[j]) <= 0.01:
                raw_touch[(pos[i], pos[j])] = True
    filled = 0
    for _r in range(rounds):
        changed = 0
        for a in range(len(geoms)):
            for b in range(a + 1, len(geoms)):
                ga, gb = geoms[a], geoms[b]
                if ga is None or gb is None:
                    continue
                d = ga.distance(gb)
                # a sliver under 1.2 m is never a street; wider needs the raw fabric's word
                limit = max_gap if raw_touch.get((a, b)) or raw_touch.get((b, a)) else 1.2
                if d <= 0.05 or d > limit:
                    continue
                tol = d + 0.5
                strip = (_topo._mitre(ga, tol).intersection(_topo._mitre(gb, tol))
                         .difference(ga).difference(gb))
                strip = strip.intersection(_uu([ga, gb]).convex_hull)
                pieces = [x for x in (strip.geoms if strip.geom_type.startswith("Multi") else [strip])
                          if not x.is_empty and x.area > 0.02]
                for piece in pieces:
                    sa = ga.buffer(0.05).intersection(piece).area
                    sb = gb.buffer(0.05).intersection(piece).area
                    k = a if sa >= sb else b
                    merged = shift_puvi._repair(_uu([geoms[k], piece]))
                    if merged is not None and not merged.is_empty:
                        geoms[k] = merged
                        changed += 1
        filled += changed
        if not changed:
            break
    gdf = gdf.copy()
    gdf["geometry"] = geoms
    return gdf, filled


def gap_pairs(gdf, raw_v, pairs, tol=0.05):
    """Neighbour pairs of the raw fabric that the layer holds apart, and the widest gap."""
    raw_geom = {i: shift_puvi._valid(g) for i, g in zip(raw_v.index, raw_v.geometry)}
    new_of = {}
    k = 0
    for i, _j in pairs:
        new_of[i] = gdf.geometry.iloc[k]
        k += 1
    idx = list(new_of)
    opened, worst = 0, 0.0
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            i, j = idx[a], idx[b]
            if raw_geom[i] is None or raw_geom[j] is None:
                continue
            if raw_geom[i].distance(raw_geom[j]) > 0.01:
                continue                      # not neighbours in the raw fabric
            d = new_of[i].distance(new_of[j])
            if d > tol:
                opened += 1
                worst = max(worst, d)
    return opened, worst


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    ap = argparse.ArgumentParser(prog="shift_refabric", description=__doc__.splitlines()[0])
    ap.add_argument("--fine", default=str(paths.PROJECT / "Puvi vector fine tunned" / "Puvi_Vector_Finetunned.geojson"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=K)
    ap.add_argument("--smooth", type=float, default=SMOOTH_M)
    args = ap.parse_args(argv)

    raw = gpd.read_file(paths.BUFFER_GPKG, layer=paths.BUFFER_LAYER).to_crs(UTM)
    fine = gpd.read_file(args.fine).to_crs(UTM)
    out_parts = []
    for v in sorted(set(fine["village_code"].astype(str))):
        raw_v = raw[raw["village_code"].astype(str).str.strip() == v]
        fine_v = fine[fine["village_code"].astype(str) == v]
        if not len(raw_v):
            out_parts.append(fine_v.copy())
            continue
        g = refabric_village(raw_v, fine_v, k=args.k, smooth=args.smooth)
        g, folds = shift_puvi.clip_siblings(g)
        pairs = pair_rows(raw_v, fine_v)
        g, filled = close_raw_neighbour_gaps(g, raw_v, pairs)
        if filled:
            g, _more = shift_puvi.clip_siblings(g)    # a filled sliver may brush a third parcel
        opened, worst = gap_pairs(g, raw_v, pairs)
        drift = [float(np.linalg.norm(np.array(g.geometry.iloc[i].centroid.coords[0])
                                      - np.array(fine_v.geometry[j].centroid.coords[0])))
                 for i, (_r, j) in enumerate(pairs)]
        log.info("%s: %d parcel(s), %d fold(s) clipped, %d gap pair(s) left (worst %.2f m), "
                 "drift from the fine positions median %.2f m max %.2f m",
                 v, len(g), folds, opened, worst,
                 float(np.median(drift)) if drift else 0.0, max(drift, default=0.0))
        out_parts.append(g)

    out = pd.concat(out_parts, ignore_index=True)
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=UTM)
    hand = {str(v) for v in out[out["hand_placed"] == 1]["village_code"]} if "hand_placed" in out.columns else set()
    layers = {v: out[out["village_code"].astype(str) == v].copy()
              for v in sorted(set(out["village_code"].astype(str)))}
    layers, notes = shift_puvi.clip_village_overlaps(layers, {v: v in hand for v in layers})
    for n in notes:
        log.info("cross-village: %s", n)
    out = gpd.GeoDataFrame(pd.concat(layers.values(), ignore_index=True), geometry="geometry", crs=UTM)
    out["area_sqm"] = out.geometry.area.round(1)
    out4326 = out.to_crs(4326)
    out4326["geometry"] = [shift_puvi._snap_or_keep(g, grid=1e-8) for g in out4326.geometry]
    out4326["geometry"] = [g if g is None or g.is_valid else g.buffer(0) for g in out4326.geometry]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out4326.to_file(args.out, driver="GeoJSON", COORDINATE_PRECISION=8)
    print("\n%d parcels written to %s" % (len(out4326), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
