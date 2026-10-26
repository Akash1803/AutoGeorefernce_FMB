import sys; sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import geopandas as gpd, shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from autogeoref import printed_numbers as pn, paths
V="35_04_221"
B = gpd.read_file("out/final_plan.gpkg", layer="plan").to_crs(32644)
table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
P = pn.printed_sets(table); AG = pn.agreed_sets(table)
S = {s: unary_union(list(B.geometry[B.survey == s])) for s in sorted(set(B.survey))}
for s in S: print(s, "prints", sorted(P.get(s, set())), "| two-way", sorted(AG.get(s, set())))
for a, b in [("265","67B"),("267","68B"),("267","70A"),("70A","73"),("70A","74"),("70A","72")]:
    print(a, b, "distance %.3f m, shared %.1f m" % (S[a].distance(S[b]), S[a].boundary.intersection(S[b].buffer(0.05)).length))
U = unary_union(list(B.geometry))
for p in getattr(U, "geoms", [U]):
    for r in p.interiors:
        h = Polygon(r).difference(U)
        if h.area > 0.01:
            near = sorted({B.survey.iat[k] + "/" + B.plot_no.iat[k] for k in range(len(B)) if B.geometry.iat[k].distance(h) < 0.05})
            mrr = h.minimum_rotated_rectangle; import numpy as np
            c = np.array(mrr.exterior.coords); e = [np.hypot(*(c[i+1]-c[i])) for i in range(2)]
            print("gap %.1f m2 at %s, mrr %.1f x %.1f m, touching %s" % (h.area, h.representative_point().wkt, max(e), min(e), near))
print(gpd.GeoSeries([S["70A"]], crs=32644).area.iat[0])
