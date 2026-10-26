from qgis.core import QgsProject, QgsGeometry
from qgis.PyQt.QtWidgets import QDockWidget, QTableView, QTableWidget
import collections
mw = iface.mainWindow()
found = []
for d in mw.findChildren(QDockWidget):
    if "topolog" in d.windowTitle().lower() or "check" in d.objectName().lower():
        found.append((d.objectName(), d.windowTitle()))
        for tv in d.findChildren(QTableView) + d.findChildren(QTableWidget):
            m = tv.model()
            if m is None: continue
            n = m.rowCount(); cols = m.columnCount()
            print("table", tv.objectName(), "rows", n, "cols", cols)
            cnt = collections.Counter(); ids = collections.defaultdict(list)
            for r in range(n):
                vals = [m.data(m.index(r, c)) for c in range(cols)]
                cnt[(str(vals[0]), str(vals[1]))] += 1
                ids[(str(vals[0]), str(vals[1]))].append(str(vals[2]) if cols > 2 else "")
            print(dict(cnt))
            for k, v in ids.items(): print(k, v[:100])
print("docks:", found)
for name in ("CMP_two_Parts", "MONI_pART"):
    L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == name][0]
    print(name, "editable:", L.isEditable(), "| modified:", L.isModified(), "| count:", L.featureCount())
    geoms = [f.geometry() for f in L.getFeatures()]
    u = QgsGeometry.unaryUnion(geoms)
    holes = 0
    for p in (u.asMultiPolygon() if u.isMultipart() else [u.asPolygon()]):
        holes += max(0, len(p) - 1)
    print(name, "union holes (QGIS):", holes, "| parts:", len(u.asMultiPolygon()) if u.isMultipart() else 1)
