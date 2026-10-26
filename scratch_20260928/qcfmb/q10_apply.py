"""Write the verified table into a layer that is open in QGIS, through one edit session (copy of the file kept first).
usage: python q10_apply.py <layer name in QGIS> <file> <final table .pkl> [DRY]"""
import os, sys, json, shutil, datetime, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd
from autogeoref import qgis_bridge as qb
HERE = os.path.dirname(os.path.abspath(__file__)); name, path, tab = sys.argv[1], sys.argv[2], sys.argv[3]; dry = len(sys.argv) > 4
F = pd.read_pickle(tab); old = json.load(open(path, encoding="utf-8")); P0 = [f["properties"] for f in old["features"]]; keep0 = list(P0[0].keys())
cols = [c for c in F.columns if not c.startswith("_")]; newf = [c for c in cols if c not in keep0]
num = {"area_acre", "area_sqm", "tamilnilam_extent_sqm", "uid_tn_avail"}
val = lambda c, v: (None if v == "" or v is None or (isinstance(v, float) and np.isnan(v)) else (float(v) if c != "uid_tn_avail" else int(v))) if c in num else ("" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v))
rows, add, delete = [], [], []
for i, g in F.groupby("_row", sort=True):
    i = int(i)
    if g._new.any():
        delete.append({"i": i, "kide": P0[i]["kide"]})
        for _j, r in g.iterrows(): add.append({"wkb": r["_geom"].wkb_hex, "attrs": {c: val(c, r[c]) for c in cols}})
    else:
        r = g.iloc[0]; a = {c: val(c, r[c]) for c in cols}; ch = {c: v for c, v in a.items() if c in newf or str(P0[i].get(c) if P0[i].get(c) is not None else "") != str(v if v is not None else "")}
        ch = {c: v for c, v in ch.items() if c in newf or c in ("kide", "survey_no", "subdiv_no", "plot_no")}          # nothing else of the existing fields is touched
        rows.append({"i": i, "kide": P0[i]["kide"], "attrs": ch})
pay = {"layer": name, "count": len(P0), "new_fields": [[c, "double" if c in num - {"uid_tn_avail"} else ("int" if c == "uid_tn_avail" else "string")] for c in newf], "rows": rows, "delete": delete, "add": add}
pf = os.path.join(HERE, "out", "apply_%s.json" % name); json.dump(pay, open(pf, "w", encoding="utf-8"))
print(name, "| rows", len(P0), "| attribute rows", len(rows), "| rows with a change to an existing field", sum(1 for r in rows if any(c not in newf for c in r["attrs"])), "| new fields", len(newf), "| rows replaced", len(delete), "by", len(add), "| after:", len(P0) - len(delete) + len(add))
if dry: sys.exit()
shutil.copy2(path, os.path.join(r"D:\Projects\Tambaram_Chengalpattu\_logs", "%s_before_%s.geojson" % (name, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))))
code = r'''
import json
from qgis.core import QgsProject, QgsGeometry, QgsFeature, QgsField, QgsWkbTypes
from qgis.PyQt.QtCore import QVariant
p = json.load(open(PAY, encoding="utf-8")); l = QgsProject.instance().mapLayersByName(p["layer"])[0]
assert not l.isEditable() and not l.isModified(), "the layer is being edited"
feats = list(l.getFeatures()); assert len(feats) == p["count"], (len(feats), p["count"])
for r in p["rows"] + p["delete"]: assert str(feats[r["i"]]["kide"]) == r["kide"], (r["i"], feats[r["i"]]["kide"], r["kide"])
T = {"string": QVariant.String, "double": QVariant.Double, "int": QVariant.Int}
l.startEditing()
try:
    for c, t in p["new_fields"]:
        if l.fields().indexOf(c) < 0: assert l.addAttribute(QgsField(c, T[t]))
    l.updateFields(); idx = {f.name(): i for i, f in enumerate(l.fields())}
    for r in p["rows"]:
        if r["attrs"]: assert l.changeAttributeValues(feats[r["i"]].id(), {idx[c]: v for c, v in r["attrs"].items()})
    assert l.deleteFeatures([feats[r["i"]].id() for r in p["delete"]]) or not p["delete"]
    new = []
    for a in p["add"]:
        g = QgsGeometry(); g.fromWkb(bytes.fromhex(a["wkb"])); assert not g.isNull() and g.isGeosValid(), "a new piece is not valid"
        if QgsWkbTypes.isMultiType(l.wkbType()) and not g.isMultipart(): g.convertToMultiType()
        if not QgsWkbTypes.isMultiType(l.wkbType()) and g.isMultipart(): assert g.convertToSingleType(), "a multi-part piece in a single-part layer"
        f = QgsFeature(l.fields()); f.setGeometry(g)
        for c, v in a["attrs"].items(): f.setAttribute(idx[c], v)
        new.append(f)
    assert l.addFeatures(new) or not new
    ok = l.commitChanges()
    if not ok: err = l.commitErrors()[:5]; l.rollBack()
    else: err = []
except Exception:
    l.rollBack(); raise
l.triggerRepaint(); print(json.dumps({"committed": ok, "errors": err, "features": l.featureCount(), "fields": len(l.fields())}))
'''.replace("PAY", repr(pf.replace("\\", "/")))
print(qb.json_result(code, timeout=900))
