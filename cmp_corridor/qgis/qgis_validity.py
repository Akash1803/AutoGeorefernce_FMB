from qgis.core import QgsProject, QgsGeometry
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == LNAME][0]
bad = []
for f in L.getFeatures():
    errs = f.geometry().validateGeometry(QgsGeometry.ValidatorGeos)
    if errs: bad.append((f["kide"], errs[0].what()[:60]))
print("QGIS (GEOS validator) invalid:", len(bad), bad[:10], "| features", L.featureCount(), "| edits pending", L.isModified())
