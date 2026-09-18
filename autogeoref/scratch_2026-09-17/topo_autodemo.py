import json, math, itertools
from qgis.core import QgsProject, QgsVectorLayer, QgsGeometry, QgsFeature, QgsField, QgsWkbTypes
from qgis.PyQt.QtCore import QVariant
import processing
import json as _json; _cfg = _json.load(open("C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/be3d1523-acf3-4fef-9ca8-783aeaecb2d3/scratchpad/georef_village.json")); VILLAGE = _cfg["village"]; SUFFIX = _cfg.get("suffix", "AUTODEMO"); RAILP = set(_cfg.get("rail_parcels", []))
p = QgsProject.instance(); grp = p.layerTreeRoot().findGroup("Georef review - " + SUFFIX + " " + VILLAGE)
layers = {n.layer().name().split(" ")[1]: n.layer() for n in grp.findLayers()}          # "AUTO 47A - ..." -> "47A"
ACRE = 0.000247105381
def check(G):
    ov = []; near = []
    for ka, kb in itertools.combinations(list(G), 2):
        if ka[0] == kb[0]: continue
        inter = G[ka].intersection(G[kb]); a = inter.area() if inter and not inter.isEmpty() else 0.0
        if a > 0.02: ov.append(["%s/%s" % ka, "%s/%s" % kb, round(a, 2)])
        d = G[ka].distance(G[kb])
        if 0 < d < 1.0 and not G[ka].intersects(G[kb]): near.append(["%s/%s" % ka, "%s/%s" % kb, round(d, 3)])
    u = QgsGeometry.unaryUnion(list(G.values()))
    if QgsWkbTypes.flatType(u.wkbType()) == QgsWkbTypes.GeometryCollection: u.convertGeometryCollectionToSubclass(QgsWkbTypes.PolygonGeometry)
    pp = u.asMultiPolygon() if u.isMultipart() else [u.asPolygon()]
    gaps = [round(QgsGeometry.fromPolygonXY([r]).area(), 2) for q in pp for r in q[1:] if QgsGeometry.fromPolygonXY([r]).area() > 0.01]
    return {"overlaps": ov, "gaps": gaps, "near_misses": near, "invalid": ["%s/%s" % k for k, g in G.items() if g.validateGeometry(QgsGeometry.ValidatorGeos)]}
G0 = {(n, f.id()): QgsGeometry(f.geometry()) for n, l in layers.items() for f in l.getFeatures()}
before = check(G0)
# keep the FMB (sheet) values as their own attributes before any geometry change
for n, l in layers.items():
    F = l.fields(); add = [QgsField(c, QVariant.Double) for c in ("fmb_area_sqm", "fmb_area_acre", "fmb_area_cent", "fmb_perimeter_m") if F.indexOf(c) < 0]
    if add: l.dataProvider().addAttributes(add); l.updateFields()
    l.startEditing(); F = l.fields()
    for f in l.getFeatures():
        a = f.geometry().area()
        for c, v in (("fmb_area_sqm", round(a, 3)), ("fmb_area_acre", round(a * ACRE, 5)), ("fmb_area_cent", round(a * ACRE * 100, 3)), ("fmb_perimeter_m", round(f.geometry().length(), 3))):
            l.changeAttributeValue(f.id(), F.indexOf(c), v)
    l.commitChanges()
# merged memory layer -> snap to itself (0.30 m) -> clip overlaps (larger yields) -> fill enclosed gaps
M = QgsVectorLayer("Polygon?crs=" + list(layers.values())[0].crs().authid(), "m", "memory")
M.dataProvider().addAttributes([QgsField("src_layer", QVariant.String), QgsField("src_fid", QVariant.Int)]); M.updateFields()
fs = []
for n, l in layers.items():
    for f in l.getFeatures():
        g = QgsGeometry(f.geometry())
        if g.validateGeometry(QgsGeometry.ValidatorGeos): g = g.makeValid()
        if QgsWkbTypes.flatType(g.wkbType()) == QgsWkbTypes.GeometryCollection: g.convertGeometryCollectionToSubclass(QgsWkbTypes.PolygonGeometry)
        if g.isMultipart() and len(g.asGeometryCollection()) > 1: g = max(g.asGeometryCollection(), key=lambda q: q.area())
        if g.isMultipart(): g.convertToSingleType()
        nf = QgsFeature(M.fields()); nf.setGeometry(g); nf.setAttributes([n, f.id()]); fs.append(nf)
M.dataProvider().addFeatures(fs)
from qgis.core import QgsProcessingContext, QgsFeatureRequest
ctx = QgsProcessingContext(); ctx.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
res = processing.run("native:snapgeometries", {"INPUT": M, "REFERENCE_LAYER": M, "TOLERANCE": 0.30, "BEHAVIOR": 0, "OUTPUT": "memory:"}, context=ctx)
G = {}
for f in res["OUTPUT"].getFeatures():
    g = QgsGeometry(f.geometry())
    if g.validateGeometry(QgsGeometry.ValidatorGeos): g = g.makeValid()
    if QgsWkbTypes.flatType(g.wkbType()) == QgsWkbTypes.GeometryCollection: g.convertGeometryCollectionToSubclass(QgsWkbTypes.PolygonGeometry)
    G[(f["src_layer"], f["src_fid"])] = g
clips = []; unresolved = []; seen_conf = set()
for _ in range(5):
    n_clip = 0
    for ka, kb in itertools.combinations(list(G), 2):
        if ka[0] == kb[0]: continue
        inter = G[ka].intersection(G[kb])
        if inter.isEmpty() or inter.area() <= 0.02: continue
        ra, rb = str(ka[0]) in RAILP, str(kb[0]) in RAILP
        if ra != rb:                                                         # railway land vs parcel: not clipped, reported for the team's decision
            k = tuple(sorted(["%s/%s" % ka, "%s/%s" % kb]))
            if k not in seen_conf: seen_conf.add(k); unresolved.append([k[0], k[1], round(inter.area(), 2)])
            continue
        big, small = (ka, kb) if G[ka].area() >= G[kb].area() else (kb, ka)
        newg = G[big].difference(G[small])
        if QgsWkbTypes.flatType(newg.wkbType()) == QgsWkbTypes.GeometryCollection: newg.convertGeometryCollectionToSubclass(QgsWkbTypes.PolygonGeometry)
        if newg.isMultipart():
            parts = sorted(newg.asGeometryCollection(), key=lambda q: -q.area()); newg = parts[0]
            for fr in parts[1:]: G[small] = G[small].combine(fr)
        G[big] = newg; n_clip += 1; clips.append(["%s/%s" % big, "clipped by", "%s/%s" % small, round(inter.area(), 2)])
    if n_clip == 0: break
u = QgsGeometry.unaryUnion(list(G.values()))
if QgsWkbTypes.flatType(u.wkbType()) == QgsWkbTypes.GeometryCollection: u.convertGeometryCollectionToSubclass(QgsWkbTypes.PolygonGeometry)
pp = u.asMultiPolygon() if u.isMultipart() else [u.asPolygon()]; filled = []
for q in pp:
    for ring in q[1:]:
        gap = QgsGeometry.fromPolygonXY([ring])
        if gap.area() <= 0.01: continue
        best, score = None, 0.0
        for k, g in G.items():
            s = g.buffer(0.02, 4).intersection(gap).area()
            if s > score: best, score = k, s
        if best: G[best] = G[best].combine(gap); filled.append([round(gap.area(), 2), "into", "%s/%s" % best])
for k in list(G):
    if G[k].validateGeometry(QgsGeometry.ValidatorGeos): G[k] = G[k].makeValid()
    if G[k].isMultipart() and len(G[k].asGeometryCollection()) > 1: G[k] = max(G[k].asGeometryCollection(), key=lambda q: q.area())
    if G[k].isMultipart(): G[k].convertToSingleType()
# write back with geometry-based measures recomputed (fmb_* fields keep the sheet values)
moved = {}
for n, l in layers.items():
    F = l.fields(); ix = {c: F.indexOf(c) for c in ("area_sqm", "area_are", "area_hect", "area_acre", "area_cent", "perimeter_m")}
    l.startEditing()
    for f in l.getFeatures():
        g = G[(n, f.id())]
        if g.equals(f.geometry()): continue
        moved["%s/%s" % (n, f.id())] = round(g.hausdorffDistance(f.geometry()), 3)
        l.changeGeometry(f.id(), g); a = g.area()
        for c, v in (("area_sqm", round(a, 3)), ("area_are", round(a / 100, 4)), ("area_hect", round(a / 1e4, 6)), ("area_acre", round(a * ACRE, 5)), ("area_cent", round(a * ACRE * 100, 3)), ("perimeter_m", round(g.length(), 3))):
            if ix[c] >= 0: l.changeAttributeValue(f.id(), ix[c], v)
    l.commitChanges(); l.triggerRepaint()
G2 = {(n, f.id()): f.geometry() for n, l in layers.items() for f in l.getFeatures()}
after = check(G2)
areas = {}
for n, l in layers.items():
    fs = list(l.getFeatures()); areas[n] = {"fmb_sqm": round(sum(f["fmb_area_sqm"] for f in fs), 2), "fmb_acre": round(sum(f["fmb_area_acre"] for f in fs), 4), "after_sqm": round(sum(f.geometry().area() for f in fs), 2), "after_acre": round(sum(f.geometry().area() for f in fs) * ACRE, 4)}
print(json.dumps({"before": before, "clips": clips, "unresolved_rail_conflicts": unresolved, "filled": filled, "moved_m": moved, "after": after, "areas": areas}))
