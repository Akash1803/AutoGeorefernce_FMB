"""Re-run the strict Tamilnilam match on a layer as it stands on disk and write only what changed (through QGIS).
usage: python q13_rematch.py <layer name in QGIS> <file> [DRY]"""
import os, sys, json, shutil, datetime, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd
import q9_tamilnilam as q9
from autogeoref import qgis_bridge as qb
HERE = os.path.dirname(os.path.abspath(__file__)); name, path = sys.argv[1], sys.argv[2]; dry = len(sys.argv) > 3
B = gpd.read_file(path, ignore_geometry=True).fillna(""); tn = [c for c in B.columns if c.startswith("tamilnilam_") or c == "uid_tn_avail"]
N = q9.match(B.drop(columns=tn)); s = lambda x: x.fillna("").astype(str).str.replace(r"\.0$", "", regex=True).replace({"nan": "", "None": ""})
ch = pd.Series(False, index=B.index)
for c in tn: ch |= s(B[c]) != s(N[c])
print(name, "| plots", len(B), "| matched before", int((s(B.uid_tn_avail) == "1").sum()), "| now", int((N.uid_tn_avail == 1).sum()), "| plots whose Tamilnilam fields change:", int(ch.sum()))
print(pd.DataFrame({"kide": B.kide, "village": B.village, "before": B.tamilnilam_match, "now": N.tamilnilam_match})[ch].head(30).to_string())
if dry or not ch.any(): sys.exit()
val = lambda c, v: (None if v == "" or v is None or (isinstance(v, float) and np.isnan(v)) else (int(v) if c == "uid_tn_avail" else float(v))) if c in ("tamilnilam_extent_sqm", "uid_tn_avail") else ("" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
pay = {"layer": name, "count": len(B), "rows": [{"i": int(i), "kide": B.kide[i], "uid": B.unique_id[i], "area": float(B.area_sqm[i]), "attrs": {c: val(c, N.at[i, c]) for c in tn}} for i in B.index[ch]]}
pf = os.path.join(HERE, "out", "rematch_%s.json" % name); json.dump(pay, open(pf, "w", encoding="utf-8"))
shutil.copy2(path, os.path.join(r"D:\Projects\Tambaram_Chengalpattu\_logs", "%s_before_%s.geojson" % (name, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))))
code = r'''
import json
from qgis.core import QgsProject
p = json.load(open(PAY, encoding="utf-8")); l = QgsProject.instance().mapLayersByName(p["layer"])[0]
assert not l.isEditable() and not l.isModified(), "the layer is being edited"; assert l.featureCount() == p["count"]
idx = {f.name(): i for i, f in enumerate(l.fields())}; feats = list(l.getFeatures())
for r in p["rows"]: assert str(feats[r["i"]]["unique_id"]) == r["uid"] and str(feats[r["i"]]["kide"]) == r["kide"] and abs(float(feats[r["i"]]["area_sqm"]) - r["area"]) < 0.05, ("row order differs", r["i"])
l.startEditing()
try:
    n = 0
    for r in p["rows"]:
        assert l.changeAttributeValues(feats[r["i"]].id(), {idx[c]: v for c, v in r["attrs"].items()}), "refused"; n += 1
    ok = l.commitChanges(); err = [] if ok else l.commitErrors()[:5]
    if not ok: l.rollBack()
except Exception:
    l.rollBack(); raise
print(json.dumps({"committed": ok, "errors": err, "features_changed": n, "features": l.featureCount()}))
'''.replace("PAY", repr(pf.replace("\\", "/")))
print(qb.json_result(code, timeout=900))
