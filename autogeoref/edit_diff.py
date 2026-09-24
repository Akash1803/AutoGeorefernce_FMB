"""Which parcels changed between two versions of a parcel file - and were they meant to?

    python -m autogeoref.edit_diff <before> <after> [--expect 54,88A] [--csv out.csv]

Akash, 2026-09-24: when one parcel's vertex was edited, other parcels (or other parts of the
same parcel) sometimes changed too, and the automated topology check that ran on the manually
georeferenced file did not show it. This compares every plot of two versions (matched by
village, survey and plot number, in file order for repeats) and reports, for each changed plot:

  * vertices moved / added / removed and the largest move;
  * whether each moved vertex was a corner it shares with another plot (a shared corner moving
    drags the neighbour with it) or its own;
  * area change;
  * whether the survey was one you expected to change (--expect); anything else is flagged.
"""
import argparse
import csv
from collections import defaultdict

import geopandas as gpd
import numpy as np
import shapely
from scipy.spatial import cKDTree


def _coords(g):
    return shapely.get_coordinates(g) if g is not None and not g.is_empty else np.zeros((0, 2))


def _key_frame(gdf):
    gdf = gdf.copy()
    if "village_code" not in gdf.columns:
        gdf["village_code"] = ""
    if "path" in gdf.columns:
        fp = gdf["path"].astype(str).str.replace("\\", "/").str.extract(r"FMB_Vector/(\d+_\d+_\d+)/")[0]
        gdf["village_code"] = gdf["village_code"].fillna(fp)
    seen = defaultdict(int)
    keys = []
    for v, s, p in zip(gdf["village_code"].astype(str), gdf["survey_no"].astype(str), gdf.get("plot_no", [""] * len(gdf))):
        base = (v, s, "" if str(p) in ("nan", "None") else str(p))
        keys.append(base + (seen[base],))
        seen[base] += 1
    gdf["_key"] = keys
    return gdf


def diff(before, after, expect=(), tol=0.005):
    """[dict] one row per plot that changed, missing or new."""
    b = _key_frame(before.to_crs(32644))
    a = _key_frame(after.to_crs(32644))
    shared = defaultdict(int)                       # how many plots carry each vertex, before
    for g in b.geometry:
        for x, y in {(round(x, 3), round(y, 3)) for x, y in _coords(g)}:
            shared[(x, y)] += 1
    bmap = dict(zip(b["_key"], b.geometry))
    amap = dict(zip(a["_key"], a.geometry))
    rows = []
    for k in sorted(set(bmap) | set(amap), key=str):
        g0, g1 = bmap.get(k), amap.get(k)
        row = {"village": k[0], "survey": k[1], "plot": k[2], "expected": k[1] in expect}
        if g0 is None or g1 is None:
            row.update(change="added" if g0 is None else "removed")
            rows.append(row)
            continue
        if g0.equals_exact(g1, tol):
            continue
        c0, c1 = _coords(g0), _coords(g1)
        d01 = cKDTree(c1).query(c0)[0] if len(c1) else np.full(len(c0), np.inf)
        d10 = cKDTree(c0).query(c1)[0] if len(c0) else np.full(len(c1), np.inf)
        moved = [(tuple(c), d) for c, d in zip(c0, d01) if d > tol]
        own = [m for m in moved if shared[(round(m[0][0], 3), round(m[0][1], 3))] <= 1]
        row.update(change="shape",
                   vertices_moved_or_removed=len({m[0] for m in moved}),
                   of_which_shared_corners=len({m[0] for m in moved}) - len({m[0] for m in own}),
                   vertices_added=int((d10 > tol).sum()),
                   largest_move_m=round(float(max((d for _c, d in moved), default=0.0)), 2),
                   area_change_m2=round(g1.area - g0.area, 2),
                   area_change_pct=round(100 * (g1.area - g0.area) / g0.area, 1) if g0.area else None)
        rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--expect", default="", help="comma-separated surveys you meant to edit")
    ap.add_argument("--csv")
    a = ap.parse_args(argv)
    expect = {s.strip() for s in a.expect.split(",") if s.strip()}
    rows = diff(gpd.read_file(a.before), gpd.read_file(a.after), expect)
    unexpected = [r for r in rows if not r["expected"]]
    print("%d plots changed; %d in surveys you did not list" % (len(rows), len(unexpected)))
    for r in rows:
        print("  %s %s/%s %s moved %s (shared corners %s) added %s max %s m area %s m2 %s%%%s" % (
            r["village"], r["survey"], r["plot"], r["change"], r.get("vertices_moved_or_removed", ""),
            r.get("of_which_shared_corners", ""), r.get("vertices_added", ""), r.get("largest_move_m", ""),
            r.get("area_change_m2", ""), r.get("area_change_pct", ""), "" if r["expected"] else "  <- not expected"))
    if a.csv:
        keys = sorted({k for r in rows for k in r})
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
