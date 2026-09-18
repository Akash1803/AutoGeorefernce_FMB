import json
from qgis.core import QgsProject, QgsMapSettings, QgsMapRendererParallelJob, QgsRectangle
from qgis.PyQt.QtCore import QSize
from qgis.PyQt.QtGui import QColor
import json as _json
_cfg = _json.load(open("C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/be3d1523-acf3-4fef-9ca8-783aeaecb2d3/scratchpad/georef_village.json")); VILLAGE = _cfg["village"]; SUFFIX = _cfg.get("suffix", "AUTODEMO")
p = QgsProject.instance(); grp = p.layerTreeRoot().findGroup("Georef review - " + SUFFIX + " " + VILLAGE)
lyrs = [n.layer() for n in grp.findLayers()]; ext = None
for l in lyrs: ext = l.extent() if ext is None else (ext.combineExtentWith(l.extent()) or ext)
lyrs += p.mapLayersByName("Google Satellite")
ms = QgsMapSettings(); ms.setDestinationCrs(p.crs()); ms.setOutputSize(QSize(1000, 900)); ms.setBackgroundColor(QColor(255, 255, 255)); ms.setLayers(lyrs)
ext.grow(20); ms.setExtent(ext); job = QgsMapRendererParallelJob(ms); job.start(); job.waitForFinished()
path = "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/be3d1523-acf3-4fef-9ca8-783aeaecb2d3/scratchpad/demo_topo_" + SUFFIX + "_" + VILLAGE + ".png"; job.renderedImage().save(path, "PNG"); print(json.dumps({"png": path}))
