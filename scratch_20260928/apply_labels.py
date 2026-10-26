"""Apply out/label_fix.json through QGIS edit sessions (attribute values only), after a backup of every file."""
import sys, os, json, sqlite3, datetime
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
from autogeoref import qgis_bridge as qb
plan = json.load(open(sys.argv[1]))
B = os.path.join(r"D:\Projects\Tambaram_Chengalpattu\_recovery_20260924", "backup_before_plot_numbers_" + datetime.datetime.now().strftime("%m%d_%H%M")); os.makedirs(B, exist_ok=True)
for key in plan:
    f = key.split("|")[0]; dst = os.path.join(B, os.path.relpath(f, r"D:\Projects\Tambaram_Chengalpattu\FMB_Vector").replace(os.sep, "__"))
    s = sqlite3.connect("file:%s?mode=ro" % f.replace(chr(92), "/"), uri=True); t = sqlite3.connect(dst); s.backup(t); t.close(); s.close()
pp = os.path.join(B, "label_fix_applied.json").replace(chr(92), "/"); json.dump(plan, open(pp, "w"), indent=1); print("backup:", B)
code = r'''
import json, os
from qgis.core import QgsProject, NULL
d = json.load(open(PAY)); res = {}; layers = {}
for l in QgsProject.instance().mapLayers().values():
    src = l.source(); f = os.path.normcase(os.path.normpath(src.split("|")[0])); ln = src.split("layername=")[1].split("|")[0] if "layername=" in src else None
    layers[(f, ln)] = l
for key in d:
    f, ln = key.split("|"); l = layers.get((os.path.normcase(os.path.normpath(f)), ln)); assert l is not None, key
    assert not l.isEditable() and not l.isModified(), key + " is being edited"
for key, rows in d.items():
    f, ln = key.split("|"); lyr = layers[(os.path.normcase(os.path.normpath(f)), ln)]; idx = {x.name(): i for i, x in enumerate(lyr.fields())}; ch = 0; bad = []
    lyr.startEditing()
    for ft in lyr.getFeatures():
        r = rows.get(str(ft.id()))
        if r is None: continue
        old = ft["plot_no"]; old = "" if old is None or old == NULL or str(old) == "NULL" else str(old)
        if old != r["old"]: bad.append([ft.id(), old, r["old"]]); continue
        lyr.changeAttributeValue(ft.id(), idx["plot_no"], r["plot_no"] if r["plot_no"] is not None else NULL)
        for k in ("numbered", "review_needed", "label_source"):
            if k in idx: lyr.changeAttributeValue(ft.id(), idx[k], r[k])
        ch += 1
    if bad: lyr.rollBack(); res[lyr.name()] = ["NOT applied", bad]
    else: res[lyr.name()] = [ch, len(rows), lyr.commitChanges()]
    lyr.triggerRepaint()
print(json.dumps(res))
'''.replace("PAY", repr(pp))
print(qb.json_result(code, timeout=300))
