from qgis.core import QgsProject
for l in QgsProject.instance().mapLayers().values():
    if l.name() in ("CMP_two_Parts", "MONI_pART") and l.isEditable():
        print(l.name(), "editable, modified:", l.isModified(), "-> commit:", l.commitChanges() if l.isModified() else l.rollBack())
print("done")
