import sys; sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
from shapely.ops import unary_union
from autogeoref import satcheck, ring_io
x0, y0, x1, y1 = 389670, 1407820, 389790, 1408040
jpg = "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/ag/out/sat_junction.jpg"
satcheck.render_from_qgis((x0, y0, x1, y1), jpg)
im = np.array(Image.open(jpg))
fig, axs = plt.subplots(1, 2, figsize=(18, 16))
L = {l["name"]: l["source"] for l in ring_io.village_layers("Paranur")}
P = gpd.read_file("out/junction_probe.gpkg").to_crs(32644)
for ax, overlay in zip(axs, (False, True)):
    ax.imshow(im, extent=(x0, x1, y0, y1))
    if overlay:
        for n in ["69_parcels_modified", "267_parcels_modified", "68B_parcels_modified", "73_parcels_auto", "266_parcels_modified"]:
            for g in ring_io._read_source(L[n]).to_crs(32644).geometry: ax.plot(*g.exterior.xy, color="white", lw=1)
        for n in ["72_parcels_auto", "70A_parcels_auto"]:
            for g in ring_io._read_source(L[n]).to_crs(32644).geometry: ax.plot(*g.exterior.xy, color="orange", lw=1, ls="--")
        for g in P.geometry: ax.plot(*g.exterior.xy, color="red", lw=1.2)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
plt.savefig("out/sat_junction.png", dpi=80, bbox_inches="tight")
