import sys, json
from qgis.core import QgsProject, QgsGeometry, QgsPointXY, QgsFeature
fx = json.load(open(FIXFILE))
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == LNAME][0]
if L.isEditable() and L.isModified(): print("uncommitted edits - stop"); raise SystemExit
F = L.fields(); ik, iv = F.indexOf("kide"), F.indexOf("village_code")
want = {(r["old_kide"], r["vc"]) for r in fx["change"] + fx["delete"]}
cand = {}
for f in L.getFeatures():
    key = (str(f[ik]), str(f[iv]))
    if key in want: cand.setdefault(key, []).append((f.id(), f.geometry()))
def find(r):
    pt = QgsGeometry.fromPointXY(QgsPointXY(*r["old_pt"])); c = [x[0] for x in cand.get((r["old_kide"], r["vc"]), []) if x[1].intersects(pt)]
    return c[0] if len(c) == 1 else None
L.startEditing(); bad = []; n = 0
for r in fx["change"]:
    fid = find(r)
    if fid is None: bad.append(r["old_kide"]); continue
    L.changeGeometry(fid, QgsGeometry.fromWkt(r["wkt"]))
    for a, v in r["attrs"].items():
        if F.indexOf(a) >= 0 and v is not None: L.changeAttributeValue(fid, F.indexOf(a), v)
    if r.get("review"): L.changeAttributeValue(fid, F.indexOf("review"), r["review"])
    n += 1
for r in fx["delete"]:
    fid = find(r)
    if fid is None: bad.append("del " + r["old_kide"]); continue
    L.deleteFeature(fid); n += 1
for a in fx["add"]:
    f = QgsFeature(F); f.setGeometry(QgsGeometry.fromWkt(a["wkt"]))
    for k, v in a["attrs"].items():
        if F.indexOf(k) >= 0 and v is not None: f.setAttribute(F.indexOf(k), v)
    L.addFeature(f); n += 1
print("commit", L.commitChanges(), "edits", n, "unmatched", bad)
