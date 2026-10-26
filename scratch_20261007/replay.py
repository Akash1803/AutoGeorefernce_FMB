"""Replay: hide some of Akash's Paranur parcels, place them from the rest, compare with his placement. Read-only."""
import sys, os, time
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.ops import unary_union
from autogeoref import ring, printed_numbers, paths
V = "35_04_221"; PR = r"D:\Projects\Tambaram_Chengalpattu"
HIS = ["66A", "66B", "67A", "67B", "68A", "68B", "69", "70A", "70B", "265", "266", "267"]
hide = sys.argv[1].split(",")
his = {s: unary_union(gpd.read_file(os.path.join(PR, s + "_parcels_modified.geojson")).to_crs(32644).geometry.values).buffer(0) for s in HIS}
seeds = {s: g for s, g in his.items() if s not in hide}
table = printed_numbers.village_table(V, os.path.join(os.path.dirname(__file__), "out", "printed_%s.json" % V))
t0 = time.time()
decisions, sheets_, ctx = ring.place_ring(V, seeds, table, size=len(hide), ring=hide)
print("time %.0f s" % (time.time() - t0))
for d in decisions:
    if d.status != "placed":
        print("%-5s ASK  %s | %s" % (d.survey, d.reason, d.metrics)); continue
    g = ring.ground_polygon(sheets_[d.survey], d.theta, d.t)
    h = his[d.survey]
    dc = g.centroid.distance(h.centroid)
    hd = shapely.hausdorff_distance(g.exterior, ring.largest(h).exterior)
    iou = g.intersection(h).area / g.union(h).area
    print("%-5s placed: centroid %.1f m from his, outline max %.1f m, IoU %.3f | %s" % (d.survey, dc, hd, iou, d.metrics["pose_from"]))
rows = []
for d in decisions:
    if d.status == "placed":
        rows.append({"survey": d.survey, "kind": "tool", "geometry": ring.ground_polygon(sheets_[d.survey], d.theta, d.t)})
for s in hide:
    rows.append({"survey": s, "kind": "his", "geometry": ring.largest(his[s])})
out = gpd.GeoDataFrame(rows, crs=32644).to_crs(4326)
out.to_file(os.path.join(os.path.dirname(__file__), "out", "replay_%s_%d.geojson" % ("_".join(hide), int(time.time()))), driver="GeoJSON")
print("written replay geojson")
