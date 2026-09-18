"""Writing GeoPackages that QGIS may be holding open, and superseding files instead of deleting.

Measured on this machine with QGIS 3.40.10: while a QgsVectorLayer holds a GeoPackage, os.replace
and os.remove both raise PermissionError(13), even when the layer is not in edit mode and not in
the layer tree. The handle is released as soon as the layer object is destroyed.
"""
import json
import os
import shutil
from pathlib import Path

from . import paths, qgis_bridge

_FIND = r"""
import json
from qgis.core import QgsProject, QgsProviderRegistry
target = %r.replace("\\", "/").lower()
proj = QgsProject.instance(); root = proj.layerTreeRoot(); out = []
for lyr in list(proj.mapLayers().values()):
    try:
        src = QgsProviderRegistry.instance().decodeUri(lyr.providerType(), lyr.source()).get("path", "")
    except Exception:
        continue
    if str(src).replace("\\", "/").lower() != target:
        continue
    node = root.findLayer(lyr.id()); parent = node.parent() if node else None
    out.append({"id": lyr.id(), "name": lyr.name(), "editable": lyr.isEditable(),
                "modified": lyr.isModified(), "uri": lyr.source(),
                "group": parent.name() if parent is not None and hasattr(parent, "name") else "",
                "index": parent.children().index(node) if parent is not None and node else -1})
print(json.dumps(out))
"""

_REMOVE = """
import gc, json
from qgis.core import QgsProject
proj = QgsProject.instance()
removed = []
for lid in %r:
    lyr = proj.mapLayer(lid)
    if lyr is not None:
        removed.append(lid); proj.removeMapLayer(lid)
gc.collect()
print(json.dumps(removed))
"""

_RESTORE = """
import json
from qgis.core import QgsProject, QgsVectorLayer
proj = QgsProject.instance(); root = proj.layerTreeRoot(); added = []
for spec in %s:
    lyr = QgsVectorLayer(spec["uri"], spec["name"], "ogr")
    if not lyr.isValid():
        continue
    proj.addMapLayer(lyr, False)
    parent = root
    if spec.get("group"):
        found = root.findGroup(spec["group"])
        parent = found if found is not None else root
    parent.insertLayer(max(0, spec.get("index", 0)), lyr)
    added.append(spec["name"])
print(json.dumps(added))
"""


def held_layers(path):
    """Every project layer backed by this file, in any group. Empty when QGIS is unreachable."""
    if not qgis_bridge.available():
        return []
    try:
        return qgis_bridge.json_result(_FIND % str(Path(path).resolve()))
    except Exception:
        return []


def release(path):
    """Remove those layers so the file can be replaced; returns what to restore afterwards."""
    held = held_layers(path)
    if not held:
        return []
    editing = [h["name"] for h in held if h["editable"]]
    if editing:
        raise PermissionError("save or stop editing %s in QGIS first" % ", ".join(editing))
    qgis_bridge.json_result(_REMOVE % [h["id"] for h in held])
    return held


def restore(path, layers):
    if not layers or not qgis_bridge.available():
        return
    try:
        qgis_bridge.json_result(_RESTORE % json.dumps(layers))
    except Exception:
        pass


def safe_write(path, writer):
    """writer(tmp_path) must create a complete file; returns False and keeps the original on failure."""
    path = Path(path)
    # keep the real extension last: OGR picks its driver from it and warns otherwise
    tmp = path.with_name(path.stem + ".autogeoref-tmp" + path.suffix)
    try:
        held = release(path)
    except PermissionError:
        return False
    try:
        if tmp.exists():
            tmp.unlink()
        writer(tmp)
        for side in ("-wal", "-shm"):
            p = Path(str(path) + side)
            if p.exists():
                p.unlink()
        os.replace(tmp, path)
        return True
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        return False
    finally:
        restore(path, held)


def supersede(path, run_id):
    """Move a file the tool is replacing into _logs; nothing is ever deleted."""
    path = Path(path)
    if not path.exists():
        return None
    dest = paths.logs_dir() / ("georef_superseded_%s" % run_id)
    dest.mkdir(parents=True, exist_ok=True)
    try:
        release(path)
    except PermissionError:
        return None
    target = dest / path.name
    shutil.move(str(path), str(target))
    return target
