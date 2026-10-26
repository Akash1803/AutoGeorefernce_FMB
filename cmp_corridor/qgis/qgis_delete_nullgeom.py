from qgis.core import QgsProject
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == LNAME][0]
if L.isEditable() and L.isModified(): print("uncommitted edits - stop"); raise SystemExit
ids = [f.id() for f in L.getFeatures() if f.geometry().isNull() or f.geometry().isEmpty()]
print(LNAME, "geometry-less features:", len(ids), [(f["kide"], f["area_sqm"]) for f in L.getFeatures(ids)] if ids else "")
if ids:
    L.startEditing(); [L.deleteFeature(i) for i in ids]; print("commit", L.commitChanges())
