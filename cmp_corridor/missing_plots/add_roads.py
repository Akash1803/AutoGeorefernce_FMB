"""Layout roads drawn on the PDF but absent from the layer (vector conversion dropped the unnumbered faces): the warped PDF faces
(label-anchored pose, k2_warped.gpkg) that carry no number are clipped to the free space inside the sheet outline (ground no
existing parcel covers) and added as unnumbered <survey> parcels. Existing parcels are not touched. Usage: add_roads.py <vc> <sv> [apply]"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, pandas as pd
from shapely import make_valid
from shapely.ops import unary_union
from shapely.geometry import Polygon
vc, sv = sys.argv[1], sys.argv[2]
A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); A["kide"] = A.kide.astype(str); A["vc"] = A.village_code.astype(str); A["sv"] = A.kide.str.split("/").str[0]
W = gpd.read_file("k2_warped.gpkg"); W = W[(W.village == vc) & (W.survey == sv)]
sheet = make_valid(unary_union(list(W.geometry)))
near = A[A.intersects(sheet.buffer(5))]; covered = make_valid(unary_union([make_valid(g) for g in near.geometry]))
free = make_valid(sheet.difference(covered.buffer(0.01)))
roads = W[W.numbers.fillna("") == ""]
pieces = []
for g in roads.geometry:
    c = make_valid(g.intersection(free))
    for z in getattr(c, "geoms", [c]):
        if z.geom_type != "Polygon" or z.area < 5: continue
        z = Polygon(z.exterior, [r for r in z.interiors if Polygon(r).area > 0.5])
        if 2 * z.area / z.length < 1.0: continue           # hairline slivers
        pieces.append(z)
P = gpd.GeoDataFrame(geometry=pieces, crs=32644); P["area"] = P.geometry.area.round(1); P["width"] = (2 * P.geometry.area / P.geometry.length).round(1)
print(f"{vc} {sv}: PDF faces {len(W)} (unnumbered {len(roads)}) | free ground inside the sheet {free.area:.0f} m2 | road pieces to add {len(P)}, {P.area.sum():.0f} m2")
print(P.sort_values("area", ascending=False)[["area", "width"]].head(12).T.to_string())
ov = sum(z.intersection(make_valid(g)).area for z in pieces for g in near.geometry if g.intersects(z)); print("overlap with existing parcels: %.3f m2" % ov)
tmpl = A[(A.vc == vc) & (A.sv == sv)].iloc[0]
adds = []
for z in pieces:
    at = {k: (None if str(v) in ("nan", "None", "NaT") else v) for k, v in tmpl.drop(labels=["geometry", "vc", "sv"]).to_dict().items()}
    for k in list(at):
        if k.startswith("tamilnilam_") or k.startswith("reginet_"): at[k] = None
    at.update({"kide": sv, "survey_no": sv, "subdiv_no": "", "unique_id": f"{vc}_{sv}/-", "area_sqm": round(z.area, 1), "area_acre": round(z.area / 4046.856, 4),
               "source": f"georeferenced from the {sv} PDF (layout road - unnumbered face missing in the vector conversion)", "status": "red - for review", "review": "rigid placement - for review",
               "pdf_check": "unnumbered face on the PDF (layout road); clipped to the ground no parcel covered", "merge_note": "", "placed_by": "auto 2026-10-26 (PDF face, label-anchored warp)"})
    adds.append({"wkt": z.wkt, "attrs": at})
json.dump({"change": [], "add": adds, "delete": []}, open(f"fix_two_roads_{sv}.json", "w"), default=str)
P.to_file(f"roads_{sv}_preview.gpkg", driver="GPKG")
