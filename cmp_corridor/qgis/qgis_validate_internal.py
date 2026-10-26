import json
from collections import Counter
from qgis.core import QgsProject, QgsGeometry
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == LNAME][0]
cnt, per = Counter(), {}
for f in L.getFeatures():
    errs = f.geometry().validateGeometry(QgsGeometry.ValidatorQgisInternal)
    for e in errs:
        w = e.what().split(" ")[0:3]; k = " ".join(w); cnt[k] += 1
        per.setdefault(f"{f['village_code']}|{f['kide']}", []).append(e.what()[:50])
print("editing:", L.isEditable(), "modified:", L.isModified())
print("features with errors:", len(per), "| errors by type:", dict(cnt.most_common(10)))
json.dump(per, open(OUT, "w"))
