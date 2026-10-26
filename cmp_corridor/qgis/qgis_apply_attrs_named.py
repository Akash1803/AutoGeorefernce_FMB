import json
from qgis.core import QgsProject, QgsGeometry, QgsPointXY
fx = json.load(open(FIXFILE))
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == LNAME][0]
if L.isEditable() and L.isModified(): print("uncommitted edits - stop"); raise SystemExit
F = L.fields(); ik, iv = F.indexOf("kide"), F.indexOf("village_code")
want = {(r["old_kide"], r["vc"]) for r in fx["change"]}
cand = {}
for f in L.getFeatures():
    key = (str(f[ik]), str(f[iv]))
    if key in want: cand.setdefault(key, []).append((f.id(), f.geometry()))
L.startEditing(); bad = []; n = 0
for r in fx["change"]:
    pt = QgsGeometry.fromPointXY(QgsPointXY(*r["old_pt"])); c = [x[0] for x in cand.get((r["old_kide"], r["vc"]), []) if x[1].intersects(pt)]
    if len(c) != 1: bad.append(r["old_kide"]); continue
    for a, v in r["attrs"].items():
        if F.indexOf(a) >= 0: L.changeAttributeValue(c[0], F.indexOf(a), v)
    n += 1
print("commit", L.commitChanges(), "features", n, "unmatched", bad)
