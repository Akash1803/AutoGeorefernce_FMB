"""Check of the any-size reader against the two villages whose numbers were verified on the PDFs today (Manual_Georeference)."""
import os, sys, warnings, collections; warnings.filterwarnings("ignore")
import pandas as pd, geopandas as gpd
import read5
R = r"D:\Projects\Tambaram_Chengalpattu"
G = gpd.read_file(R + r"\Manual_Georeferenced\Manual_Georeference.geojson", ignore_geometry=True).fillna("")
sc = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
tot = hit = extra = 0; miss = []; ext = []
for (v, s), g in G.groupby(["village_code", "survey_no"]):
    want = collections.Counter(x for x, st in zip(g.subdiv_no, g.subdiv_status) if x and st == "printed on the FMB sheet")
    if not want: continue
    pdf = os.path.join(R, "FMB_Sketches", v, s + ".pdf")
    if not os.path.exists(pdf): continue
    L, S = read5.read(pdf, sc.get((v, s), 500)); got = collections.Counter(a["text"] for a in read5.plot_numbers(L))
    for k, n in want.items():
        tot += n; h = min(n, got.get(k, 0)); hit += h
        if h < n: miss.append((v[-3:], s, k))
    e = [k for k in got if k not in want and len(k) >= 2]; extra += len(e)
    if e: ext.append((v[-3:], s, e))
print("verified numbers", tot, "| read by the any-size reader", hit, "| missed", tot - hit); print("missed:", miss); print("numbers of 2+ characters read but not in the verified layer:", extra, ext)
