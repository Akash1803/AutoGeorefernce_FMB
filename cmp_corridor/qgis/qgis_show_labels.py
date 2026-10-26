import json
from qgis.core import QgsProject, QgsRuleBasedRenderer, QgsFillSymbol, QgsSymbol
from qgis.utils import iface
L = [l for l in QgsProject.instance().mapLayers().values() if l.name() == "CMP_two_Parts"][0]
keys = {tuple(k) for k in json.load(open(KF))}
# rule renderer: today's label/geometry changes in orange, rest as before (pink)
root = QgsRuleBasedRenderer.Rule(None)
r1 = QgsRuleBasedRenderer.Rule(QgsFillSymbol.createSimple({"color": "255,140,0,200", "outline_color": "120,40,0,255", "outline_width": "0.5"}), 0, 0,
                               "\"pdf_check\" LIKE '%2026-10-26%' OR \"merge_note\" LIKE '%2026-10-26%' OR \"placed_by\" LIKE '%2026-10-26%' OR \"pdf_check\" LIKE '%vector conversion%'", "changed today (label / missing subdivision)")
r2 = QgsRuleBasedRenderer.Rule(QgsFillSymbol.createSimple({"color": "243,166,178,120", "outline_color": "35,35,35,255", "outline_width": "0.26"}), 0, 0, "ELSE", "other parcels")
root.appendChild(r1); root.appendChild(r2)
L.setRenderer(QgsRuleBasedRenderer(root)); L.triggerRepaint()
ids = [f.id() for f in L.getFeatures() if (str(f["village_code"]), str(f["kide"])) in keys]
L.selectByIds(ids); iface.setActiveLayer(L); iface.mapCanvas().zoomToSelected(L); iface.mapCanvas().refresh()
print("highlighted (rule) + selected:", len(ids))
