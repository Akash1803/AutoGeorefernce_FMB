import sys; sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import geopandas as gpd
from autogeoref import ring_io
rows = []
for L in ring_io.village_layers("Paranur"):
    n = L["name"]; s = n.split("_parcels")[0]
    if "_parcels" not in n or (s in ("70A", "70B") and not n.endswith("_auto")): continue
    g = ring_io._read_source(L["source"]).to_crs(4326)
    for pos in range(len(g)):
        rows.append({"layer": n, "survey": s, "pos": pos, "plot_no": str(g.plot_no.iat[pos]), "geometry": g.geometry.iat[pos]})
gpd.GeoDataFrame(rows, geometry="geometry", crs=4326).to_file("out/disk_plan.gpkg", layer="plan", driver="GPKG")
