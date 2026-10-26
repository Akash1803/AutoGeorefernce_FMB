"""Final deliverable = CMP_two_Parts (live, EPSG:32644) + the fixed moni part, one GeoJSON, one schema (kide/survey_no/subdiv_no/area first),
`dataset` tells which part a row came from. The two hairline overlaps on the seam (Kattour 199 vs Sempakkam 233/234) are cut from 199.
Usage: merge_final.py <moni_fixed.geojson> <out.geojson>"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, pandas as pd, numpy as np
from shapely import make_valid
from shapely.ops import unary_union
from shapely.strtree import STRtree
MONI, OUT = sys.argv[1], sys.argv[2]
assert OUT.lower().endswith(".geojson") and "two_parts" not in OUT.lower() and "moni_part" not in OUT.lower()
A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); A["dataset"] = "CMP_two_Parts (Akash)"
M = gpd.read_file(MONI).to_crs(32644); M["dataset"] = "MONI_pART (team)"
lead = ["kide", "survey_no", "subdiv_no", "area_acre", "area_sqm", "unique_id", "district", "taluk", "village", "village_code"]
tn = [c for c in A.columns if c.startswith("tamilnilam_")]; rg = [c for c in A.columns if c.startswith("reginet_")]
tail = ["source", "status", "placed_by", "pdf_check", "merge_note", "review", "layer", "path", "poly_id", "giscode", "sheet_scale", "source_pdf", "dataset"]
cols = lead + tn + rg + tail
for c in cols:
    for X in (A, M):
        if c not in X: X[c] = None
F = gpd.GeoDataFrame(pd.concat([A[cols + ["geometry"]], M[cols + ["geometry"]]], ignore_index=True), geometry="geometry", crs="EPSG:32644")
# seam: overlaps between the two parts -> cut the Akash-side team parcel (199 blob) by the moni parcel
GA = [make_valid(g) for g in F.geometry]; T = STRtree(GA); nA = len(A); cut = 0
for i in range(nA):
    for j in T.query(GA[i], predicate="intersects"):
        if j < nA: continue
        x = GA[i].intersection(GA[j])
        if x.area > 1e-6:
            GA[i] = make_valid(GA[i].difference(GA[j].buffer(0.0005, join_style=2))); cut += 1
            GA[i] = max([z for z in getattr(GA[i], "geoms", [GA[i]]) if z.geom_type == "Polygon"], key=lambda z: z.area)
F = F.set_geometry(GA)
F["area_sqm"] = F.geometry.area.round(1); F["area_acre"] = (F.geometry.area / 4046.856).round(4)
F["survey_no"] = F.kide.astype(str).str.split("/").str[0]; F["subdiv_no"] = F.kide.astype(str).map(lambda k: k.split("/")[1] if "/" in k else "")
F.to_file(OUT, driver="GeoJSON")
print(f"rows {len(F)} = {nA} + {len(M)} | seam cuts {cut} | invalid {int((~F.is_valid).sum())} | null kide {int(F.kide.isna().sum())} | -> {OUT}")
print("villages:", F.village.value_counts().to_dict())
