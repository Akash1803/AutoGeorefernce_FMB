"""Rebuild of the last re-placement batch on the TEAM poses (offline, then one diff for QGIS):
start = mypart_v1_before_edgematch.geojson (state before edge matching), drop the rows of the given surveys, insert each rigid
sheet <vc>_<sv>_team.gpkg with replace_in_merged.py (NOFILL: it cuts team CMP parcels it covers, yields to sheets already placed),
then diff against the LIVE layer -> fix_rebuild.json. Usage: python rebuild_batch.py vc:sv vc:sv ... (order = placement order)"""
import sys, os, json, shutil, subprocess
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd
from shapely import make_valid
LIVE = r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_my_Part.geojson"
TMP = os.path.abspath("CMP_my_Part_rebuild_tmp.geojson")
B = gpd.read_file("mypart_v1_before_edgematch.geojson")
todo = [a.split(":") for a in sys.argv[1:]]
key = B.village_code.astype(str) + "|" + B.kide.astype(str).str.split("/").str[0]
drop = key.isin([f"{vc}|{sv}" for vc, sv in todo])
# team parcels around these surveys get their delivered team shape back (cuts/stretches by the misplaced sheets undone)
from shapely.ops import unary_union
T1 = gpd.read_file("merged_v1_backup.geojson"); T1["kide"] = T1.kide.astype(str); T1["vc"] = T1.village_code.astype(str)
T1m = T1.to_crs(32644); Bm = B.to_crs(32644)
foot = make_valid(unary_union([g for g, v, k in zip(T1m.geometry, T1m.vc, T1m.kide) if f"{v}|{k.split('/')[0]}" in {f"{a}|{b}" for a, b in todo}])).buffer(30)
isteam = ~B.source.astype(str).str.startswith("georef") & ~B.status.astype(str).str.startswith("green")
restored = []
for i in B.index[isteam & ~drop & Bm.intersects(foot)]:
    c = T1m[(T1m.vc == str(B.village_code[i])) & (T1m.kide == str(B.kide[i]))]
    if not len(c): continue
    j = max(c.index, key=lambda j: c.geometry[j].intersection(make_valid(Bm.geometry[i])).area)
    if not T1.geometry[j].equals_exact(B.geometry[i], 1e-9):
        B.loc[i, "geometry"] = T1.geometry[j]; restored.append(str(B.kide[i]))
        B.loc[i, "merge_note"] = "2026-10-23: team shape restored (it had been cut/stretched around a misplaced sheet)"
print("team parcels restored to their delivered shape:", len(restored), restored[:40])
B[~drop].to_file(TMP, driver="GeoJSON")
print("start rows", len(B), "| dropped", int(drop.sum()))
env = dict(os.environ, NOFILL="1", PYTHONIOENCODING="utf-8")
for vc, sv in todo:
    r = subprocess.run([sys.executable, "replace_in_merged.py", TMP.replace("\\", "/").replace("CMP_my_Part_rebuild_tmp", "CMP_my_Part_rebuild_tmp"), vc, sv, f"{vc}_{sv}_team.gpkg",
                        "auto 2026-10-23 night: rigid PDF sheet posed on the team CMP footprint of the survey (the imagery placement had moved it 7-35 m and opened gaps)"],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", env=env)
    keep = [l for l in r.stdout.splitlines() if any(k in l for k in ("cut where", "CMP parcels cut", "removed (>=", "written", "invalid"))]
    print(f"--- {vc} {sv}:", " | ".join(k[:160] for k in keep) or r.stderr[-400:])
# diff TMP vs LIVE
C = gpd.read_file(LIVE); T = gpd.read_file(TMP)
for D in (C, T): D["key"] = D.village_code.astype(str) + "|" + D.kide.astype(str)
ch, adds, dels = [], [], []
cols = ("area_sqm", "area_acre", "merge_note", "source", "placed_by", "status", "pdf_check", "subdiv_no", "unique_id", "review")
for k, g in T.groupby("key"):
    c = C[C.key == k]
    if not len(c):
        adds += [a for _, a in g.iterrows()]; continue
    used = set()
    for _, a in g.iterrows():
        cand = [(i, rr) for i, rr in c.iterrows() if i not in used]
        if not cand: adds.append(a); continue
        i, b = max(cand, key=lambda t: make_valid(t[1].geometry).intersection(make_valid(a.geometry)).area if make_valid(t[1].geometry).intersects(make_valid(a.geometry)) else -t[1].geometry.distance(a.geometry))
        used.add(i)
        if not a.geometry.equals_exact(b.geometry, 1e-9) or any(str(a.get(x)) != str(b.get(x)) for x in cols if x in T.columns):
            rp = make_valid(b.geometry).representative_point()
            ch.append({"vc": str(a.village_code), "old_kide": str(b.kide), "old_pt": [rp.x, rp.y], "wkt": a.geometry.wkt,
                       "attrs": {x: (None if a[x] is None or a[x] != a[x] else (str(a[x]) if not isinstance(a[x], (int, float)) else a[x])) for x in cols if x in T.columns},
                       "review": "rigid placement - for review" if str(a.source).startswith("georef") else None})
    for i, b in c.iterrows():
        if i not in used:
            rp = make_valid(b.geometry).representative_point(); dels.append({"vc": str(b.village_code), "old_kide": str(b.kide), "old_pt": [rp.x, rp.y]})
for k in set(C.key) - set(T.key):
    for _, b in C[C.key == k].iterrows():
        rp = make_valid(b.geometry).representative_point(); dels.append({"vc": str(b.village_code), "old_kide": str(b.kide), "old_pt": [rp.x, rp.y]})
out = {"change": ch, "add": [{"wkt": a.geometry.wkt, "attrs": {x: (None if a[x] is None or a[x] != a[x] else a[x]) for x in T.columns if x not in ("geometry", "key")}} for a in adds], "delete": dels}
for a in out["add"]: a["attrs"]["review"] = "rigid placement - for review" if str(a["attrs"].get("source", "")).startswith("georef") else a["attrs"].get("review")
json.dump(out, open("fix_rebuild.json", "w"), default=str)
print(f"rebuild rows {len(T)} vs live {len(C)} | changes {len(ch)}, adds {len(adds)}, deletes {len(dels)} -> fix_rebuild.json")
