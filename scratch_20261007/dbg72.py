"""Explore 72's candidate positions with 73/74 (red) and his parcels as the placed set."""
import sys
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from autogeoref import ring, ring_io, printed_numbers as pn, paths
seeds, src, plots = ring_io.load_seeds("Paranur", include_red=True)
table = pn.village_table("35_04_221", str(paths.PROJECT / "Auto_Georef" / "35_04_221" / "printed_numbers.json"))
ctx = ring.Context(seeds, table, seeds)
sp = ring.seed_rotations("35_04_221", seeds, table); ctx.poses.update(sp)
for k in sp: ctx.outlines[k] = ring.load_sheet("35_04_221", k, table).polygon
ctx.set_rotation(-6.6, 6.8)
sh = ring.load_sheet("35_04_221", "72", table); ctx.outlines["72"] = sh.polygon
named = [n for n in ctx.ids if ctx.named("72", n) == "agreed"]
print("named placed:", named, "| seed poses:", sorted(sp))
rows = []
for n in named:
    for th, t, src_ in ring.candidate_poses(sh, ctx, n):
        if not ctx.rotation_ok(th):
            continue
        m = ring.evaluate(sh, th, t, ctx)
        th2, t2 = ring.refine(sh, th, t, ctx, named)
        m2 = ring.evaluate(sh, th2, t2, ctx)
        rows.append((m2["valid"], m2["total"], {k: round(v) for k, v in m2["support"].items()}, round(m2["overlap"], 3),
                     m2["label_conflicts"][:2], m2["label_far"][:2], list(m2["contra"])[:3], src_, round(th2, 1)))
west = [r for r in rows if r[2].get("73", 0) + r[2].get("74", 0) > 150]
print("candidates", len(rows), "| with >150 m on 73/74:", len(west))
for r in sorted(west, key=lambda r: -(r[2].get("73", 0) + r[2].get("74", 0)))[:8]:
    print("  valid %s total %.0f %s overlap %s conflicts %s far %s contra %s | %s theta %s" % r)
print("best valid overall:")
for r in sorted([r for r in rows if r[0]], key=lambda r: -r[1])[:4]:
    print("  valid %s total %.0f %s overlap %s conflicts %s far %s contra %s | %s theta %s" % r)
