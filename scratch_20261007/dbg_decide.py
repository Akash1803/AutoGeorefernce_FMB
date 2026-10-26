"""Debug decide(): candidates, their refinement and scores for one sheet in a given placed set."""
import sys, os
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.ops import unary_union
from autogeoref import ring, printed_numbers, match, fit
V = "35_04_221"; PR = r"D:\Projects\Tambaram_Chengalpattu"
HIS = ["66A", "66B", "67A", "67B", "68A", "68B", "69", "70A", "70B", "265", "266", "267"]
seeds_ids = sys.argv[1].split(","); s = sys.argv[2]
his = {x: unary_union(gpd.read_file(os.path.join(PR, x + "_parcels_modified.geojson")).to_crs(32644).geometry.values).buffer(0) for x in HIS}
table = printed_numbers.village_table(V, os.path.join(os.path.dirname(__file__), "out", "printed_%s.json" % V))
ctx = ring.Context({x: his[x] for x in seeds_ids}, table, seeds_ids)
sh = ring.load_sheet(V, s, table)
vp = ring.clockwise(ring.largest(his[s]).exterior.coords)
ch = match.common_chains(sh.outline, vp, tol_abs=3, tol_rel=0.1, min_len=5)[0]
th0, t0, rms, mx = fit.rigid_fit([a for a, b in ch.pairs], [b for a, b in ch.pairs])
P0 = ring.pose_points(sh.outline, th0, t0)
named = [n for n in ctx.ids if ctx.named(s, n) == "agreed"]
print("named placed:", named)
rows = []
for n in named:
    for th, t, src in ring.candidate_poses(sh, ctx, n):
        m = ring.evaluate(sh, th, t, ctx)
        d = np.abs(ring.pose_points(sh.outline, th, t) - P0).max()
        rows.append((d, th, t, src, m))
rows.sort(key=lambda r: r[0])
print("candidates: %d, valid %d" % (len(rows), sum(r[4]["valid"] for r in rows)))
for d, th, t, src, m in rows[:8]:
    th2, t2 = ring.refine(sh, th, t, ctx, named); m2 = ring.evaluate(sh, th2, t2, ctx)
    d2 = np.abs(ring.pose_points(sh.outline, th2, t2) - P0).max()
    print("  corner diff %.1f m (%s) total %.0f valid %s %s| refined: diff %.1f total %.0f valid %s %s %s %s" % (
        d, src, m["total"], m["valid"], (m["label_far"], m["label_conflicts"], m["contra"], round(m["overlap"], 3)) if not m["valid"] else "",
        d2, m2["total"], m2["valid"], {k: round(v) for k, v in m2["support"].items()}, m2["label_far"], round(m2["overlap"], 3)))
dec = ring.decide(sh, ctx)
print("decide:", dec.status, dec.reason, dec.metrics)
d = np.abs(ring.pose_points(sh.outline, dec.theta, dec.t) - P0).max(); print("decision corner diff vs his rigid pose %.1f m (his own rigid fit rms %.1f)" % (d, rms))
