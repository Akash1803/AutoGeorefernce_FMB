"""Debug: why does no pose of a sheet pass (seeds = his 12 Paranur parcels)."""
import sys, os, collections
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd
from shapely.ops import unary_union
from autogeoref import ring, printed_numbers
V = "35_04_221"; PR = r"D:\Projects\Tambaram_Chengalpattu"
HIS = ["66A", "66B", "67A", "67B", "68A", "68B", "69", "70A", "70B", "265", "266", "267"]
s = sys.argv[1]
seeds = {x: unary_union(gpd.read_file(os.path.join(PR, x + "_parcels_modified.geojson")).to_crs(32644).geometry.values).buffer(0) for x in HIS}
table = printed_numbers.village_table(V, PR + r"\Auto_Georef\35_04_221\printed_numbers.json")
ctx = ring.Context(seeds, table, seeds)
ctx.set_rotation(-6.6, 6.8); print("window +-", ring.ROT_WINDOW_DEG + 6.8)
sh = ring.load_sheet(V, s, table)
print("labels:", [(L.text, round(L.x), round(L.y)) for L in sh.labels])
print("outline edges:", [round(np.hypot(*(np.subtract(sh.outline[(i + 1) % len(sh.outline)], sh.outline[i]))), 1) for i in range(len(sh.outline))])
named = [n for n in ctx.ids if ctx.named(s, n) == "agreed"]
print("named placed:", named)
why = collections.Counter(); n_rot = 0; best = []
for n in named:
    for th, t, src in ring.candidate_poses(sh, ctx, n):
        if not ctx.rotation_ok(th):
            n_rot += 1; continue
        m = ring.evaluate(sh, th, t, ctx)
        if m["valid"]:
            why["valid"] += 1
        else:
            if m["overlap"] > ring.MAX_OVERLAP: why["overlap"] += 1
            if m["label_conflicts"]: why["label: " + m["label_conflicts"][0]] += 1
            if m["label_far"]: why["far: " + m["label_far"][0].split(" ")[0]] += 1
            if m["contra"]: why["contra: " + ",".join(m["contra"])] += 1
        best.append((m["total"], th, src, round(m["overlap"], 3), m["label_conflicts"][:2], m["label_far"][:2], m["contra"]))
print("rotation-rejected:", n_rot, "| reasons:", why.most_common(12))
for b in sorted(best, key=lambda b: -b[0])[:8]:
    print("  total %.0f theta %.1f %s overlap %s conflicts %s far %s contra %s" % b)
