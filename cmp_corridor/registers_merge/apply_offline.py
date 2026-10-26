"""Apply a fix json ({change:[{vc, old_kide, old_pt, wkt?, attrs?}], add:[{wkt, attrs}], delete:[{vc, old_kide, old_pt}]}) to a GeoJSON offline
(same matching as qgis_apply_two: village_code + kide, then the feature containing old_pt). Usage: apply_offline.py <in.geojson> <fix.json> <out.geojson>"""
import sys, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, pandas as pd
from shapely import wkt as swkt
from shapely.geometry import Point
IN, FIX, OUT = sys.argv[1:4]
assert OUT != IN and "CMP_two_Parts.geojson" not in OUT and "MONI_pART" not in OUT, "output must be a new file"
P = gpd.read_file(IN); fx = json.load(open(FIX))
def find(vc, kide, pt):
    k = P.kide.astype(str).where(P.kide.notna(), "None"); c = P[(P.village_code.astype(str) == str(vc)) & (k == str(kide))]
    if len(c) == 0 and str(kide) in ("None", "nan", ""): c = P[P.kide.isna()]
    if len(c) == 1: return c.index[0]
    p = Point(pt); inside = c[c.contains(p)]
    if len(inside) == 1: return inside.index[0]
    if len(c): 
        j = c.distance(p).idxmin()
        if c.loc[j].geometry.distance(p) < 0.5: return j
    return None
n_ch = n_del = 0; miss = []
for ch in fx.get("change", []):
    i = find(ch["vc"], ch["old_kide"], ch["old_pt"])
    if i is None: miss.append(("change", ch["vc"], ch["old_kide"])); continue
    if ch.get("wkt"): P.loc[i, "geometry"] = swkt.loads(ch["wkt"])
    for k, v in (ch.get("attrs") or {}).items():
        if k not in P.columns: P[k] = None
        P.loc[i, k] = v
    n_ch += 1
dels = []
for d in fx.get("delete", []):
    i = find(d["vc"], d["old_kide"], d["old_pt"])
    if i is None: miss.append(("delete", d["vc"], d["old_kide"])); continue
    dels.append(i)
P = P.drop(index=dels); n_del = len(dels)
adds = []
for a in fx.get("add", []):
    row = dict(a["attrs"]); row["geometry"] = swkt.loads(a["wkt"]); adds.append(row)
if adds: P = pd.concat([P, gpd.GeoDataFrame(adds, geometry="geometry", crs=P.crs)], ignore_index=True)
P = P.reset_index(drop=True)
P.to_file(OUT, driver="GeoJSON")
print(f"changed {n_ch} | deleted {n_del} | added {len(adds)} | missing {miss} | rows {len(P)} -> {OUT}")
