import json
from qgis.core import QgsProject, QgsGeometry, QgsPointXY
fx = json.load(open(FIXFILE))
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == "CMP_two_Parts"][0]
F = L.fields(); cand = [(f.id(), f.geometry(), str(f["village_code"])) for f in L.getFeatures() if f["kide"] is None or str(f["kide"]) in ("", "NULL", "nan", "None")]
print("null-kide features in layer:", len(cand))
L.startEditing(); n = 0; bad = []
for r in fx["change"]:
    pt = QgsGeometry.fromPointXY(QgsPointXY(*r["old_pt"])); c = [fid for fid, g, v in cand if g.intersects(pt)]
    if len(c) != 1: bad.append(r["vc"]); continue
    L.changeGeometry(c[0], QgsGeometry.fromWkt(r["wkt"]))
    for a, v in r["attrs"].items():
        if F.indexOf(a) >= 0 and v is not None: L.changeAttributeValue(c[0], F.indexOf(a), v)
    n += 1
print("commit", L.commitChanges(), "edits", n, "unmatched", bad)
