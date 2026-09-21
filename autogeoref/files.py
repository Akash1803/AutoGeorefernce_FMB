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
from qgis.PyQt.QtXml import QDomDocument
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
    path = []
    q = parent
    while q is not None and q is not root and hasattr(q, "name"):
        path.insert(0, q.name()); q = q.parent()
    doc = QDomDocument(); lyr.exportNamedStyle(doc)
    out.append({"id": lyr.id(), "name": lyr.name(), "editable": lyr.isEditable(),
                "modified": lyr.isModified(), "uri": lyr.source(), "visible": node.itemVisibilityChecked() if node else True,
                "group_path": path, "index": parent.children().index(node) if parent is not None and node else -1,
                "style_xml": doc.toString()})
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
from qgis.PyQt.QtXml import QDomDocument
proj = QgsProject.instance(); root = proj.layerTreeRoot(); added, failed = [], []
for spec in %s:
    lyr = QgsVectorLayer(spec["uri"], spec["name"], "ogr")
    if not lyr.isValid():
        failed.append(spec["name"]); continue
    if spec.get("style_xml"):
        doc = QDomDocument(); doc.setContent(spec["style_xml"]); lyr.importNamedStyle(doc)
    proj.addMapLayer(lyr, False)
    parent = root
    for name in spec.get("group_path", []):            # the exact group, by path, never by a bare name
        nxt = next((c for c in parent.children() if c.nodeType() == 0 and c.name() == name), None)
        if nxt is None:
            break
        parent = nxt
    node = parent.insertLayer(max(0, spec.get("index", 0)), lyr)
    if node is not None:
        node.setItemVisibilityChecked(bool(spec.get("visible", True)))
    added.append(spec["name"])
print(json.dumps({"added": added, "failed": failed}))
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
    """Put released layers back where they were, with their style. Never silent: on 2026-09-21 a
    failed restore dropped ten of the team's layers from the project without a word."""
    if not layers:
        return {"added": [], "failed": []}
    if not qgis_bridge.available():
        raise RuntimeError("QGIS went away before %d layer(s) of %s could be restored: %s"
                           % (len(layers), path, ", ".join(l["name"] for l in layers)))
    out = qgis_bridge.json_result(_RESTORE % json.dumps(layers))
    if out.get("failed"):
        raise RuntimeError("could not restore layer(s) %s for %s" % (", ".join(out["failed"]), path))
    return out


def safe_write(path, writer, allow_release=False):
    """writer(tmp_path) must create a complete file; returns False and keeps the original on failure.

    While QGIS holds the file the only way to replace it is to take its layers out of the project
    and put them back. By default that is refused and the write returns False, so nobody's layer is
    ever removed from their project by accident (2026-09-21: ten of the team's layers were dropped
    this way). The engine passes allow_release=True for its own output files only.
    """
    path = Path(path)
    # keep the real extension last: OGR picks its driver from it and warns otherwise
    tmp = path.with_name(path.stem + ".autogeoref-tmp" + path.suffix)
    held = []
    try:
        if held_layers(path) and not allow_release:
            return False                           # the team's layer stays; the caller reports it
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
