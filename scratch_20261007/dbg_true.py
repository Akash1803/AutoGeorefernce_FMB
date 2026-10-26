"""Debug: evaluate the pose that matches Akash's own placement, and the candidates near it."""
import sys, os
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.ops import unary_union
from autogeoref import ring, printed_numbers, match, fit
V = "35_04_221"; PR = r"D:\Projects\Tambaram_Chengalpattu"
HIS = ["66A", "66B", "67A", "67B", "68A", "68B", "69", "70A", "70B", "265", "266", "267"]
hide = sys.argv[1].split(","); s = sys.argv[2]
his = {x: unary_union(gpd.read_file(os.path.join(PR, x + "_parcels_modified.geojson")).to_crs(32644).geometry.values).buffer(0) for x in HIS}
seeds = {x: g for x, g in his.items() if x not in hide}
table = printed_numbers.village_table(V, os.path.join(os.path.dirname(__file__), "out", "printed_%s.json" % V))
ctx = ring.Context(seeds, table, seeds)
sh = ring.load_sheet(V, s, table)
print("labels:", [(L.text, round(L.x), round(L.y)) for L in sh.labels])
vp = ring.clockwise(ring.largest(his[s]).exterior.coords)
ch = match.common_chains(sh.outline, vp, tol_abs=3, tol_rel=0.1, min_len=5)[0]
th, t, rms, mx = fit.rigid_fit([a for a, b in ch.pairs], [b for a, b in ch.pairs])
print("true-ish pose from his own outline: theta %.2f rms %.2f max %.2f (chain %d pairs)" % (th, rms, mx, ch.n_pairs))
m = ring.evaluate(sh, th, t, ctx)
print({k: v for k, v in m.items() if k != "dirs"})
cands = []
for n in [x for x in ctx.ids if ctx.named(s, x) == "agreed"]:
    for th2, t2, src in ring.candidate_poses(sh, ctx, n):
        c0 = ring.ground_polygon(sh, th2, t2).centroid; c1 = ring.ground_polygon(sh, th, t).centroid
        cands.append((c0.distance(c1), abs(((th2 - th) + 180) % 360 - 180), src, th2, t2))
cands.sort(key=lambda c: c[0] + c[1])
print("nearest candidates to the true pose:")
for d, dth, src, th2, t2 in cands[:6]:
    m2 = ring.evaluate(sh, th2, t2, ctx)
    print("  %.1f m %.1f deg %s -> total %.0f valid %s overlap %.3f contra %s conflicts %s far %s" % (d, dth, src, m2["total"], m2["valid"], m2["overlap"], m2["contra"], m2["label_conflicts"], m2["label_far"]))
print("--- refine from the true-ish pose:")
named = [x for x in ctx.ids if ctx.named(s, x) == "agreed"]
th3, t3 = ring.refine(sh, th, t, ctx, named)
m3 = ring.evaluate(sh, th3, t3, ctx)
g3 = ring.ground_polygon(sh, th3, t3); h = his[s]
print("refined: dtheta %.2f shift %.2f m | total %.0f %s valid %s overlap %.3f far %s conflicts %s contra %s | IoU vs his %.3f" % (
    th3 - th, np.linalg.norm(np.asarray(t3) - np.asarray(t)), m3["total"], {k: round(v) for k, v in m3["support"].items()}, m3["valid"], m3["overlap"], m3["label_far"], m3["label_conflicts"], m3["contra"],
    g3.intersection(h).area / g3.union(h).area))
g0 = ring.ground_polygon(sh, th, t); print("true-ish IoU vs his %.3f" % (g0.intersection(h).area / g0.union(h).area))
