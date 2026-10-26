"""Akash chose option 1 (2026-10-07): the PDF-exact 266/267/68B/69 replace his versions; his neighbours and the red
72/73/74 give way by at most 3 m so the seams close. Stage 1 (this script, read-only on his files): compute the
new geometries and write them to OUT as a plan; nothing of his is written here.
Usage: python accept_pdf4.py <out.gpkg>"""
import sys, os, json
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from autogeoref import ring_io, neat, ring_topo, printed_numbers as pn, paths
OUT = sys.argv[1]
assert os.path.abspath(OUT).lower().startswith(os.path.abspath(os.path.dirname(__file__)).lower())
ACCEPT = ["266", "267", "68B", "69"]
HOLD = ["70A", "70B"]          # need more than 3 m: not moved here, re-placed from their PDFs next
MAX_MOVE_M = 3.0
neat.GAP_PDF_M = MAX_MOVE_M                       # how far a vertex may move onto a PDF neighbour
table = pn.village_table("35_04_221", str(paths.PROJECT / "Auto_Georef" / "35_04_221" / "printed_numbers.json"))
agreed = pn.agreed_sets(table)
rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]
    if "_parcels" not in n:
        continue
    s = n.split("_parcels")[0]
    red = n.endswith("_parcels_auto")
    if s in ACCEPT and not red:
        role = "his_replaced"                     # his old version: replaced by the red PDF version
    elif s in ACCEPT and red:
        role = "accepted"                         # becomes his geometry, fixed
    elif red:
        role = "red_mover"
    else:
        role = "his_mover"
    g = ring_io._read_source(L["source"])
    for pos in range(len(g)):
        rows.append({"layer": n, "survey": s, "role": role, "pos": pos, "source": L["source"],
                     "plot_no": str(g.plot_no.iat[pos]) if "plot_no" in g else "", "geometry": g.geometry.iat[pos]})
A = gpd.GeoDataFrame(rows, geometry="geometry", crs=4326)
print(A.groupby("role").survey.unique())
W = A[A.role != "his_replaced"].reset_index(drop=True)
orig_ll = list(W.geometry)
U = W.to_crs(32644)
U["survey_no"] = U.survey
cleaned, rep = neat.clean(U, {k: set(v) for k, v in agreed.items()}, fixed_surveys=ACCEPT + HOLD)
C = cleaned.to_crs(4326)
# untouched vertices back to their exact original coordinates (the UTM round trip moves them by ~1e-9 deg)
pts = np.vstack([shapely.get_coordinates(g) for g in orig_ll])
tree = cKDTree(pts)
def back(gm):
    def f(c):
        d, k = tree.query(c)
        c2 = c.copy(); near = d <= 1e-8; c2[near] = pts[k[near]]; return c2
    return shapely.transform(gm, f)
new_ll = [orig_ll[i] if (W.role.iat[i] == "accepted" or W.survey.iat[i] in HOLD) else back(C.geometry.iat[i]) for i in range(len(W))]
# two accepted PDF sheets that overlap each other (their drawings disagree): the smaller gives way along the other's line
acc = [i for i in range(len(W)) if W.role.iat[i] == "accepted"]
for x in acc:
    for y in acc:
        if x < y and W.survey.iat[x] != W.survey.iat[y] and new_ll[x].intersects(new_ll[y]):
            inter = shapely.intersection(new_ll[x], new_ll[y])
            if sum(p.area for p in getattr(inter, "geoms", [inter]) if p.geom_type == "Polygon") > 0:
                small, big = (x, y) if new_ll[x].area < new_ll[y].area else (y, x)
                new_ll[small] = ring_topo._largest(shapely.difference(new_ll[small], new_ll[big]))
W2 = W.copy(); W2["geometry"] = new_ll
fixed = W2[(W2.role == "accepted") | W2.survey.isin(HOLD)]; movers = W2[(W2.role != "accepted") & ~W2.survey.isin(HOLD)]
movers = ring_topo.final_pass(fixed, movers.reset_index(drop=True))
W3 = gpd.GeoDataFrame(pd.concat([fixed, movers], ignore_index=True), geometry="geometry", crs=4326)
# change table
before = {(r.layer, r.pos): r.geometry for r in W.itertuples()}
out = []
for r in W3.itertuples():
    b = before[(r.layer, r.pos)]
    bu = gpd.GeoSeries([b], crs=4326).to_crs(32644).iat[0]; au = gpd.GeoSeries([r.geometry], crs=4326).to_crs(32644).iat[0]
    out.append({"layer": r.layer, "survey": r.survey, "role": r.role, "pos": r.pos, "plot_no": r.plot_no,
                "area_before": round(bu.area, 2), "area_after": round(au.area, 2), "change_%": round(100 * (au.area / bu.area - 1), 2) if bu.area else 0,
                "max_move_m": round(shapely.hausdorff_distance(bu, au), 2), "changed": not b.equals_exact(r.geometry, 0),
                "valid": r.geometry.is_valid and r.geometry.geom_type == "Polygon"})
T = pd.DataFrame(out)
pd.set_option("display.width", 220)
print(T[T.changed].to_string(index=False))
bad = T[(T.role != "accepted") & ~T.survey.isin(HOLD) & T.changed & ((T.max_move_m > MAX_MOVE_M + 0.05) | (T["change_%"].abs() > 15) | ~T.valid)]
print("\nbeyond limits:", len(bad)); print(bad.to_string(index=False) if len(bad) else "")
# check 1 on the plan: all plots in layer CRS
g = W3.geometry.values; K2 = 110540.0 * 108600.0
t = shapely.STRtree(g); a_, b_ = t.query(g, predicate="intersects"); k = a_ < b_
ov = [(W3.survey.iat[i], W3.survey.iat[j], round(sum(p.area for p in getattr(x, "geoms", [x]) if p.geom_type == "Polygon") * K2, 4))
      for i, j in zip(a_[k], b_[k]) if W3.survey.iat[i] != W3.survey.iat[j] for x in [shapely.intersection(g[i], g[j])]]
ov = [o for o in ov if o[2] > 0]
Uall = unary_union(g)
holes = [round(shapely.Polygon(r).difference(Uall).area * K2, 3) for p in getattr(Uall, "geoms", [Uall]) for r in p.interiors]
print("plan check: overlaps", len(ov), round(sum(o[2] for o in ov), 3), ov[:10], "| holes", len(holes), holes)
W3.drop(columns=[]).to_file(OUT, layer="plan", driver="GPKG")
T.to_csv(OUT.replace(".gpkg", "_changes.csv"), index=False)
print("plan written", OUT)
