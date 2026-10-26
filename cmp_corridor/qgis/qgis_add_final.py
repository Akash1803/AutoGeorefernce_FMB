import os
from qgis.core import QgsProject, QgsVectorLayer, QgsRuleBasedRenderer, QgsFillSymbol, QgsPalLayerSettings, QgsTextFormat, QgsVectorLayerSimpleLabeling
from qgis.PyQt.QtGui import QColor, QFont
prj = QgsProject.instance()
old = [l for l in prj.mapLayers().values() if l.name() == LNAME]
if old:
    L = old[0]; L.dataProvider().reloadData(); L.triggerRepaint(); print("reloaded existing layer")
else:
    L = QgsVectorLayer(PATH, LNAME, "ogr"); assert L.isValid(), "layer failed to load"; prj.addMapLayer(L); print("added layer", L.featureCount())
def sym(fill, outline, w="0.26"): return QgsFillSymbol.createSimple({"color": fill, "outline_color": outline, "outline_width": w})
root = QgsRuleBasedRenderer.Rule(None)
rules = [("\"status\" LIKE 'green%'", "approved / hand-made (green)", sym("120,200,120,150", "0,100,0,255")),
         ("\"placed_by\" LIKE '%2026-10-26%' OR \"pdf_check\" LIKE '%2026-10-26%' OR \"merge_note\" LIKE '%2026-10-26%'", "changed / added today (26 Oct) - review", sym("255,140,0,200", "120,40,0,255", "0.5")),
         ("\"status\" LIKE 'red%'", "georeferenced sheets (red - for review)", sym("240,120,120,130", "160,0,0,255")),
         ("\"dataset\" LIKE 'MONI%'", "team part (moni)", sym("170,190,240,120", "40,40,120,255")),
         ("ELSE", "other parcels", sym("243,166,178,120", "35,35,35,255"))]
for expr, label, s in rules:
    root.appendChild(QgsRuleBasedRenderer.Rule(s, 0, 0, expr, label))
L.setRenderer(QgsRuleBasedRenderer(root))
pal = QgsPalLayerSettings(); pal.fieldName = "kide"; pal.enabled = True
tf = QgsTextFormat(); tf.setFont(QFont("Arial", 7)); tf.setSize(7); tf.setColor(QColor(20, 20, 20)); pal.setFormat(tf)
pal.scaleVisibility = True; pal.minimumScale = 4000; pal.maximumScale = 0
L.setLabelsEnabled(True); L.setLabeling(QgsVectorLayerSimpleLabeling(pal)); L.triggerRepaint()
iface.setActiveLayer(L); iface.mapCanvas().setExtent(L.extent()); iface.mapCanvas().refresh()
print("styled + zoomed:", L.name(), L.featureCount(), "features, crs", L.crs().authid())
