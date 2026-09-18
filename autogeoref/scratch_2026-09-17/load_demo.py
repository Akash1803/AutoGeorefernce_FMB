import json, os, glob
from qgis.core import (QgsProject, QgsVectorLayer, QgsFillSymbol, QgsSingleSymbolRenderer, QgsLayerTreeGroup, QgsPalLayerSettings, QgsVectorLayerSimpleLabeling, QgsTextFormat,
                       QgsMapSettings, QgsMapRendererParallelJob, QgsRectangle)
from qgis.PyQt.QtCore import QSize
from qgis.PyQt.QtGui import QColor
p = QgsProject.instance(); root = p.layerTreeRoot()
import json as _json
_cfg = _json.load(open("C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/be3d1523-acf3-4fef-9ca8-783aeaecb2d3/scratchpad/georef_village.json")); VILLAGE = _cfg["village"]; SUFFIX = _cfg.get("suffix", "AUTODEMO")
D = "D:/Projects/Tambaram_Chengalpattu/FMB_Vector/" + VILLAGE
GROUP = "Georef review - " + SUFFIX + " " + VILLAGE
old = root.findGroup(GROUP)
if old:
    for n in old.findLayers(): p.removeMapLayer(n.layerId())
    root.removeChildNode(old)
grp = QgsLayerTreeGroup(GROUP); root.insertChildNode(0, grp)
added = []; ext = None
for f in sorted(glob.glob(os.path.join(D, "*_parcels_" + SUFFIX + ".gpkg"))):
    sv = os.path.basename(f).split("_")[0]
    l = QgsVectorLayer(f + "|layername=parcels", "AUTO %s" % sv, "ogr")
    conf = next(l.getFeatures())["georef_confidence"]; method = next(l.getFeatures())["georef_method"]
    col = "#2ecc40" if conf >= 80 else ("#ffb000" if conf >= 50 else "#ff4136")
    if method == "seed": col = "#0074d9"
    l.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple({"color": col + "66", "outline_color": col, "outline_width": "0.8"})))
    s = QgsPalLayerSettings(); s.fieldName = "concat('%s/', \"plot_no\", ' (', \"georef_confidence\", ')')" % sv; s.isExpression = True; s.enabled = True
    fmt = QgsTextFormat(); fmt.setSize(9); fmt.setColor(QColor("black")); s.setFormat(fmt); l.setLabeling(QgsVectorLayerSimpleLabeling(s)); l.setLabelsEnabled(True)
    l.setName("AUTO %s - %s, confidence %s" % (sv, method, conf))
    p.addMapLayer(l, False); grp.addLayer(l); added.append((sv, method, conf))
    e = l.extent(); ext = e if ext is None else (ext.combineExtentWith(e) or ext)
# team's hand-placed versions as dashed black outlines for comparison
hand_names = [l.name() for l in p.mapLayers().values() if any(l.name() == sv or l.name().startswith(sv + "_parcels_modified") for sv, _, _ in added)]
for n in hand_names:
    for l in p.mapLayersByName(n):
        l.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple({"color": "0,0,0,0", "outline_color": "#000000", "outline_width": "0.5", "outline_style": "dash"})))
        l.setLabelsEnabled(False); l.triggerRepaint()
sat = p.mapLayersByName("Google Satellite"); lyrs = [n.layer() for n in grp.findLayers()] + [p.mapLayersByName(n)[0] for n in hand_names if p.mapLayersByName(n)] + sat
ms = QgsMapSettings(); ms.setDestinationCrs(p.crs()); ms.setOutputSize(QSize(1000, 900)); ms.setBackgroundColor(QColor(255, 255, 255)); ms.setLayers(lyrs)
ext.grow(25); ms.setExtent(ext); job = QgsMapRendererParallelJob(ms); job.start(); job.waitForFinished()
path = "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/be3d1523-acf3-4fef-9ca8-783aeaecb2d3/scratchpad/demo_review_" + SUFFIX + "_" + VILLAGE + ".png"; job.renderedImage().save(path, "PNG")
print(json.dumps({"added": added, "png": path}))
