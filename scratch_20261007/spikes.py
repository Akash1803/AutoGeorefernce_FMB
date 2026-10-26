"""Read-only: find spikes, hooks and tiny zigzags in every Paranur parcel QGIS shows (his + red)."""
import sys, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from autogeoref import ring_io
rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]
    if "_parcels" not in n:
        continue
    g = ring_io._read_source(L["source"])
    U = g.to_crs(32644)
    for pos in range(len(g)):
        geom = U.geometry.iat[pos]
        if geom is None or geom.is_empty:
            continue
        for poly in getattr(geom, "geoms", [geom]):
            c = np.array(poly.exterior.coords)[:-1]
            m = len(c)
            for i in range(m):
                a, b, d = c[i - 1], c[i], c[(i + 1) % m]
                e1, e2 = np.linalg.norm(b - a), np.linalg.norm(d - b)
                if e1 < 1e-9 or e2 < 1e-9:
                    rows.append((n, pos, i, "repeated corner", 0.0, 0.0, b)); continue
                cosang = np.dot(a - b, d - b) / (e1 * e2)
                ang = math.degrees(math.acos(max(-1, min(1, cosang))))       # interior angle at b
                # spike: a sharp turn back on a short edge
                if ang < 25 and min(e1, e2) < 6:
                    rows.append((n, pos, i, "spike", round(ang, 1), round(min(e1, e2), 2), b))
                # hook / zigzag: a short edge between two near-reversals
                elif min(e1, e2) < 1.0 and ang < 120:
                    rows.append((n, pos, i, "small jog", round(ang, 1), round(min(e1, e2), 2), b))
D = pd.DataFrame(rows, columns=["layer", "pos", "vertex", "kind", "angle_deg", "short_edge_m", "xy"])
pd.set_option("display.width", 200)
print(D.groupby(["layer", "kind"]).size().unstack(fill_value=0).to_string())
print(D.drop(columns="xy").to_string(index=False))
D["x"] = [p[0] for p in D.xy]; D["y"] = [p[1] for p in D.xy]
D.drop(columns="xy").to_csv(r"C:\Users\FAI-AK~1\AppData\Local\Temp\claude\d--Akash-Python\93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9\scratchpad\ag\out\spikes.csv", index=False)
