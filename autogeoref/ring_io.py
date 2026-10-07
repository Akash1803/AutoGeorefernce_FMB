"""Inputs and outputs of a ring run: the seeds as QGIS shows them, the batch files, the QGIS group.

Seeds are the layers in Manual_Georeferenced/<Village>_Village of Akash's open QGIS project (read
over the QGIS bridge, so they are exactly what he sees). Layers named *_parcels_auto are the tool's
own unapproved output and are never seeds. Nothing of his is written here.
"""
import datetime
import json
import os

import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

from . import paths, qgis_bridge

GROUP = "Manual_Georeferenced"
RED = {"style": "no", "outline_color": "255,32,32,255", "outline_width": "0.7"}


def group_name(village_name):
    return "%s_Village" % village_name


def village_layers(village_name):
    """[{name, source, editing}] in Manual_Georeferenced/<Village>_Village of the open project."""
    code = """
import json
from qgis.core import QgsProject
root = QgsProject.instance().layerTreeRoot()
g = root.findGroup(%r)
sub = g.findGroup(%r) if g else None
out = []
if sub:
    for n in sub.findLayers():
        L = n.layer()
        if L is not None:
            out.append({"name": L.name(), "source": L.source(), "editing": L.isEditable()})
print(json.dumps(out))
""" % (GROUP, group_name(village_name))
    return qgis_bridge.json_result(code)


def _read_source(source):
    path, _, rest = source.partition("|")
    layer = None
    for part in rest.split("|"):
        if part.startswith("layername="):
            layer = part.split("=", 1)[1]
    return gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)


def load_seeds(village_name, crs=32644, include_red=False):
    """({survey: dissolved polygon in UTM}, {survey: source}, plots GeoDataFrame in EPSG:4326 with survey_no)
    for every hand / approved layer of the village group."""
    seeds, sources, plots = {}, {}, []
    for L in village_layers(village_name):
        name = L["name"]
        if "_parcels" not in name or (name.endswith("_parcels_auto") and not include_red):
            continue
        survey = name.split("_parcels")[0]
        g = _read_source(L["source"])
        g = g[g.geometry.notna() & ~g.geometry.is_empty]
        if g.empty:
            continue
        plots.append(gpd.GeoDataFrame({"survey_no": [survey] * len(g)}, geometry=g.to_crs(4326).geometry.values, crs=4326))
        seeds[survey] = unary_union(g.to_crs(crs).geometry.buffer(0).values)
        sources[survey] = L["source"]
    plots = gpd.GeoDataFrame(pd.concat(plots, ignore_index=True), geometry="geometry", crs=4326) if plots else None
    return seeds, sources, plots


def batch_dir(village):
    base = paths.PROJECT / "Auto_Georef" / village
    base.mkdir(parents=True, exist_ok=True)
    used = [int(p.name.split("_")[1]) for p in base.glob("batch_*") if p.is_dir() and p.name.split("_")[1].isdigit()]
    n = 1 + max(used, default=0)
    d = base / ("batch_%02d_%s" % (n, datetime.datetime.now().strftime("%Y%m%d_%H%M")))
    d.mkdir()
    return d


def output_path(village, survey):
    return paths.vector_dir(village) / ("%s_parcels_auto.gpkg" % survey)


def write_survey(village, survey, gdf_utm):
    """The placed plots as FMB_Vector/<v>/<s>_parcels_auto.gpkg, layer 'parcels', EPSG:4326."""
    out = output_path(village, survey)
    out_gdf = gdf_utm if gdf_utm.crs and gdf_utm.crs.to_epsg() == 4326 else gdf_utm.to_crs(4326)   # no round trip
    if out.exists():
        os.remove(out)
    out_gdf.to_file(out, layer="parcels", driver="GPKG")
    return out


def qgis_show(village_name, files):
    """Load each file as <s>_parcels_auto (red outline) at the top of the village group; an earlier
    unapproved layer of the same survey is replaced. Approved layers are never touched."""
    code = """
import json
from qgis.core import QgsProject, QgsVectorLayer, QgsFillSymbol, QgsSingleSymbolRenderer
P = QgsProject.instance(); root = P.layerTreeRoot()
g = root.findGroup(%(group)r) or root.addGroup(%(group)r)
sub = g.findGroup(%(sub)r) or g.addGroup(%(sub)r)
done = []
for survey, path in %(files)r:
    name = survey + "_parcels_auto"
    for L in P.mapLayersByName(name):
        P.removeMapLayer(L.id())
    L = QgsVectorLayer(path + "|layername=parcels", name, "ogr")
    if not L.isValid():
        done.append([survey, "invalid"]); continue
    red = dict(%(red)r)
    hand = [n.layer() for n in sub.findLayers() if n.layer() is not None and not n.layer().name().endswith("_auto")
            and "_parcels" in n.layer().name()]
    if hand:                                       # same line width as his own layers (rule of 2026-09-20)
        sl = hand[0].renderer().symbol().symbolLayer(0).properties()
        red["outline_width"] = sl.get("outline_width", red["outline_width"])
        red["outline_width_unit"] = sl.get("outline_width_unit", "MM")
    L.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple(red)))
    P.addMapLayer(L, False); sub.insertLayer(0, L)
    done.append([survey, L.featureCount()])
print(json.dumps(done))
""" % {"group": GROUP, "sub": group_name(village_name), "files": [(s, str(p).replace("\\", "/")) for s, p in files],
       "red": RED}
    return qgis_bridge.json_result(code)


def qgis_show_options(village_name, gpkg, layer_name):
    """The asked sheets' candidate positions, dashed yellow and labelled '<s> option <k>', at the top of the group."""
    code = """
import json
from qgis.core import (QgsProject, QgsVectorLayer, QgsFillSymbol, QgsSingleSymbolRenderer, QgsPalLayerSettings,
                       QgsTextFormat, QgsTextBufferSettings, QgsVectorLayerSimpleLabeling)
from qgis.PyQt.QtGui import QColor
P = QgsProject.instance(); root = P.layerTreeRoot()
g = root.findGroup(%(group)r) or root.addGroup(%(group)r)
sub = g.findGroup(%(sub)r) or g.addGroup(%(sub)r)
for L in P.mapLayersByName(%(name)r):
    P.removeMapLayer(L.id())
L = QgsVectorLayer(%(path)r + "|layername=ask_options", %(name)r, "ogr")
L.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple(
    {"style": "no", "outline_color": "255,220,0,255", "outline_width": "0.6", "outline_style": "dash"})))
s = QgsPalLayerSettings(); s.fieldName = "label"; tf = QgsTextFormat(); tf.setSize(9); tf.setColor(QColor(255, 220, 0))
b = QgsTextBufferSettings(); b.setEnabled(True); b.setSize(0.8); b.setColor(QColor(0, 0, 0)); tf.setBuffer(b); s.setFormat(tf)
L.setLabeling(QgsVectorLayerSimpleLabeling(s)); L.setLabelsEnabled(True)
P.addMapLayer(L, False); sub.insertLayer(0, L)
print(json.dumps([L.name(), L.featureCount()]))
""" % {"group": GROUP, "sub": group_name(village_name), "name": layer_name, "path": str(gpkg).replace(os.sep, "/")}
    return qgis_bridge.json_result(code)


def clear_unapproved(village, surveys):
    """Drop an earlier run's red <s>_parcels_auto layers (and their files) for these surveys before a re-run.
    A file shown under any other layer name (approved) is left alone and reported."""
    code = """
import json
from qgis.core import QgsProject
P = QgsProject.instance(); out = {"removed": [], "kept": []}
for s in %r:
    for L in P.mapLayers().values():
        src = L.source().split("|")[0].replace("\\\\", "/").lower()
        if src.endswith("/%%s_parcels_auto.gpkg" %% s.lower()) and "/fmb_vector/%s/" in src:
            if L.name() == s + "_parcels_auto":
                out["removed"].append(L.name()); P.removeMapLayer(L.id())
            else:
                out["kept"].append(L.name())
import gc, os
from qgis.core import QgsApplication
gc.collect(); QgsApplication.processEvents()          # QGIS keeps the GeoPackage open until then
kept = {n.split("_parcels")[0] for n in out["kept"]}
import time
for s, p in %r:
    for attempt in range(10):                       # the OGR connection pool lets go after a moment
        if s in kept or not os.path.exists(p):
            break
        try:
            os.remove(p)
        except PermissionError:
            gc.collect(); QgsApplication.processEvents(); time.sleep(0.5)
print(json.dumps(out))
""" % (list(surveys), village.lower(), [(s, str(output_path(village, s))) for s in surveys])
    return qgis_bridge.json_result(code)


def write_report(d, rows):
    df = pd.DataFrame(rows)
    df.to_csv(d / "batch_report.csv", index=False)
    try:
        df.to_excel(d / "batch_report.xlsx", index=False)
    except Exception:
        pass
    return df


def save_run(d, info):
    (d / "run.json").write_text(json.dumps(info, indent=1, default=str), encoding="utf-8")
