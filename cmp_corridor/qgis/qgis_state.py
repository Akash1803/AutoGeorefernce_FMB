from qgis.core import QgsProject
for l in QgsProject.instance().mapLayers().values():
    if l.name() in ("CMP_two_Parts", "CMP_Final_Deliverable", "MONI_pART"):
        print(l.name(), "| editable:", l.isEditable(), "| modified:", l.isModified(), "| features:", l.featureCount())
