"""Apply out/close_plan.json {file|layer: {fid: wkt}} through QGIS edit sessions (the layers are open there):
backup of every file first (sqlite backup, includes unflushed edits), geometry + area/perimeter/vertex fields,
a note, and the plot's rows in the file's "edges" layer when the file has one."""
import sys, os, json, math, sqlite3, datetime, shapely, pyproj
from shapely.geometry import Polygon
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
from autogeoref import qgis_bridge as qb
NOTE = sys.argv[2] if len(sys.argv) > 2 else "seams closed 2026-09-28"
plan = json.load(open(sys.argv[1])); proj = pyproj.Proj("EPSG:32644")
B = os.path.join(r"D:\Projects\Tambaram_Chengalpattu\_recovery_20260924", "backup_before_village_close_" + datetime.datetime.now().strftime("%m%d_%H%M")); os.makedirs(B, exist_ok=True)
pay = {}
for key, d in plan.items():
    f, layer = key.split("|"); dst = os.path.join(B, os.path.basename(f))
    if not os.path.exists(dst):
        s = sqlite3.connect("file:%s?mode=ro" % f.replace(chr(92), "/"), uri=True); t = sqlite3.connect(dst); s.backup(t); t.close(); s.close()
    pay[key] = {}
    for fid, w in d.items():
        ring = list(shapely.from_wkt(w).exterior.coords)
        if Polygon(ring).exterior.is_ccw: ring = ring[::-1]
        ring = ring[:-1]; st = max(range(len(ring)), key=lambda k: ring[k][1]); ring = ring[st:] + ring[:st]; E = []
        for k in range(len(ring)):
            (x0, y0), (x1, y1) = ring[k], ring[(k + 1) % len(ring)]; grid = math.degrees(math.atan2(x1 - x0, y1 - y0)) % 360
            lon, lat = proj((x0 + x1) / 2, (y0 + y1) / 2, inverse=True); conv = proj.get_factors(lon, lat).meridian_convergence
            E.append({"edge_no": k + 1, "length_m": round(math.dist((x0, y0), (x1, y1)), 3), "bearing_grid_deg": round(grid, 2), "bearing_true_deg": round((grid + conv) % 360, 2),
                      "grid_convergence_deg": round(conv, 3), "wkt": shapely.LineString([(x0, y0), (x1, y1)]).wkt})
        pay[key][fid] = {"wkt": w, "edges": E}
pp = os.path.join(B, "plan_applied.json"); json.dump(pay, open(pp, "w")); print("backup:", B, "| files", len(pay))
code = r'''
import json, os
from qgis.core import QgsProject, QgsGeometry, QgsVectorLayer, QgsFeature
d = json.load(open(PAY)); res = {}; layers = {}
for l in QgsProject.instance().mapLayers().values():
    src = l.source(); f = os.path.normcase(os.path.normpath(src.split("|")[0])); ln = src.split("layername=")[1].split("|")[0] if "layername=" in src else None
    layers[(f, ln)] = l
for key in d:
    f, ln = key.split("|"); l = layers.get((os.path.normcase(os.path.normpath(f)), ln))
    assert l is not None, key
    assert not l.isEditable() and not l.isModified(), key + " is being edited"
for key, plots in d.items():
    f, ln = key.split("|"); lyr = layers[(os.path.normcase(os.path.normpath(f)), ln)]
    idx = {x.name(): i for i, x in enumerate(lyr.fields())}; lyr.startEditing(); ch = 0; done = {}
    for ft in lyr.getFeatures():
        p = plots.get(str(ft.id()))
        if p is None: continue
        geom = QgsGeometry.fromWkt(p["wkt"])
        if ft.geometry().isMultipart(): geom.convertToMultiType()
        lyr.changeGeometry(ft.id(), geom); a = geom.area(); ch += 1
        done[str(ft.id())] = (ft["poly_id"] if "poly_id" in idx else None, ft["plot_no"] if "plot_no" in idx else None, ft["survey_no"] if "survey_no" in idx else None)
        for k, v in {"area_sqm": round(a, 3), "area_are": round(a / 100, 4), "area_hect": round(a / 10000, 6), "area_acre": round(a / 4046.8564224, 5),
                     "area_cent": round(a / 40.468564224, 3), "perimeter_m": round(geom.length(), 3), "n_vertices": sum(1 for _ in geom.vertices()) - 1}.items():
            if k in idx: lyr.changeAttributeValue(ft.id(), idx[k], v)
        if "georef_notes" in idx and NOTE not in str(ft["georef_notes"] or ""):
            old = ft["georef_notes"]; old = "" if old is None or str(old) == "NULL" else str(old)
            lyr.changeAttributeValue(ft.id(), idx["georef_notes"], (old + "; " if old else "") + NOTE)
    ok = lyr.commitChanges(); res[lyr.name()] = [ch, len(plots), ok]
    sub = [s.split("!!::!!")[1] for s in lyr.dataProvider().subLayers()] if lyr.dataProvider().subLayers() else []
    if False and all(v[0] is not None for v in done.values()):
        le = QgsVectorLayer(f + "|layername=edges", "tmp_edges", "ogr"); le.startEditing(); flds = le.fields(); pids = {v[0] for v in done.values()}
        le.deleteFeatures([x.id() for x in le.getFeatures() if x["poly_id"] in pids]); feats = []
        for fid, (pid, pn, sv) in done.items():
            for e in plots[fid]["edges"]:
                nf = QgsFeature(flds); nf.setGeometry(QgsGeometry.fromWkt(e["wkt"]))
                for k, v in {"survey_no": sv, "poly_id": pid, "plot_no": pn, **{kk: e[kk] for kk in ("edge_no", "length_m", "bearing_grid_deg", "bearing_true_deg", "grid_convergence_deg")}}.items():
                    if flds.indexOf(k) >= 0: nf.setAttribute(k, v)
                feats.append(nf)
        le.addFeatures(feats); res[lyr.name()].append(["edges", len(feats), le.commitChanges()]); del le
    lyr.triggerRepaint()
print(json.dumps(res))
'''.replace("PAY", repr(pp.replace(chr(92), "/"))).replace("NOTE", repr(NOTE))
r = qb.json_result(code, timeout=300); print(r)
bad = {k: v for k, v in r.items() if v[0] != v[1] or not v[2]}; print("problems:", bad)
