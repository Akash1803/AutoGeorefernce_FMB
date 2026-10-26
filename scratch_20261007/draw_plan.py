"""Draw a plan (survey outlines, plot numbers, gaps in yellow) to a PNG. Usage: draw_plan.py <plan.gpkg> <out.png> [survey ...]"""
import sys, warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import unary_union
B = gpd.read_file(sys.argv[1], layer="plan").to_crs(32644)
focus = sys.argv[3:]
U = unary_union(list(B.geometry))
gaps = [Polygon(r).difference(U) for p in getattr(U, "geoms", [U]) for r in p.interiors]
fig, ax = plt.subplots(figsize=(14, 14)); ax.set_aspect("equal")
for k in range(len(B)):
    g = B.geometry.iat[k]; red = B.layer.iat[k].endswith("_auto")
    x, y = g.exterior.xy; ax.plot(x, y, color="red" if red else "black", lw=1.2 if red else 0.6)
    p = g.representative_point(); ax.text(p.x, p.y, B.plot_no.iat[k], fontsize=7, ha="center", color="red" if red else "black", clip_on=True)
for s in sorted(set(B.survey)):
    d = unary_union(list(B.geometry[B.survey == s])); x, y = d.exterior.xy if d.geom_type == "Polygon" else d.convex_hull.exterior.xy
    ax.plot(x, y, color="blue", lw=2, alpha=0.5); c = d.centroid; ax.text(c.x, c.y, s, fontsize=14, color="blue", ha="center", alpha=0.6, clip_on=True)
for h in gaps:
    for q in getattr(h, "geoms", [h]):
        if q.geom_type == "Polygon": ax.fill(*q.exterior.xy, color="yellow", alpha=0.9, ec="orange")
if focus:
    d = unary_union(list(B.geometry[B.survey.isin(focus)])).buffer(15); x0, y0, x1, y1 = d.bounds; ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect("equal", adjustable="datalim")
plt.savefig(sys.argv[2], dpi=90, bbox_inches="tight")
