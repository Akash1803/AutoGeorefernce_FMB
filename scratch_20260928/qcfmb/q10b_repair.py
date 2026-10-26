"""Repair after the first write: 9 pieces had several parts and the layer (single polygons) kept only the first part.
The missing parts are added as their own features - same ground as in the original layer, nothing moved."""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from autogeoref import qgis_bridge as qb
HERE = os.path.dirname(os.path.abspath(__file__)); LAYER = r"D:/Projects/Tambaram_Chengalpattu/FMB_Georeferenced_Final/FMB_Fine_Tunned.geojson"
S = pd.read_csv(os.path.join(HERE, "out", "split_pieces_v1.csv"), dtype=str).fillna(""); S["g"] = [shapely.from_wkb(bytes.fromhex(w)) for w in S.wkb]; S = S[[len(shapely.get_parts(g)) > 1 for g in S.g]]
B = gpd.read_file(LAYER); wb = {g.wkb: i for i, g in enumerate(B.geometry.normalize())}; sqm = lambda g: float(gpd.GeoSeries([g], crs=4326).to_crs(32644).area.iloc[0])
fix, add = [], []
for r in S.itertuples():
    parts = list(shapely.get_parts(r.g)); hit = [k for k, p in enumerate(parts) if shapely.normalize(p).wkb in wb]; assert len(hit) == 1, (r.old_kide, r.pid, hit)
    i = wb[shapely.normalize(parts[hit[0]]).wkb]; a = B.iloc[i]; a0 = sqm(parts[hit[0]]); big = [p for k, p in enumerate(parts) if k != hit[0] and sqm(p) >= 5.0]; n = 1 + len(big)
    att = {c: (None if (isinstance(a[c], float) and np.isnan(a[c])) or a[c] is None or a[c] is pd.NaT else (a[c].isoformat() if hasattr(a[c], "isoformat") else (a[c].item() if hasattr(a[c], "item") else a[c]))) for c in B.columns if c != "geometry"}
    m = att["tamilnilam_match"]; m = (m.split(" (")[0] + " (subdivision drawn in %d pieces)" % n) if m.startswith("matched") and n > 1 else m
    fix.append({"kide": att["kide"], "area_old": float(att["area_sqm"]), "wkb": shapely.normalize(parts[hit[0]]).wkb_hex, "attrs": {"area_sqm": round(a0, 1), "area_acre": round(a0 / 4046.8564224, 4), "tamilnilam_match": m}})
    for k, p in enumerate(parts):
        if k == hit[0]: continue
        ap = sqm(p); d = dict(att); d.update({"area_sqm": round(ap, 1), "area_acre": round(ap / 4046.8564224, 4), "tamilnilam_match": m})
        if ap < 5.0:
            d.update({"kide": att["survey_no"], "subdiv_no": "", "plot_no": "", "unique_id": att["village_code"] + "_" + att["survey_no"] + "/-", "subdiv_status": "sliver left by the cut of the merged row (belongs with %s)" % att["kide"], "tamilnilam_match": "not appended: the plot has no subdivision number on the sheet", "uid_tn_avail": 0})
            for c in d:
                if c.startswith("tamilnilam_") and c != "tamilnilam_match": d[c] = None if c == "tamilnilam_extent_sqm" else ""
        add.append({"wkb": p.wkb_hex, "attrs": d}); print("  part added:", d["kide"], round(ap, 1), "m2 |", d["subdiv_status"][:60])
print("features to correct", len(fix), "| parts to add", len(add), "| ground given back: %.1f m2" % sum(x["attrs"]["area_sqm"] for x in add))
pf = os.path.join(HERE, "out", "repair.json"); json.dump({"fix": fix, "add": add}, open(pf, "w", encoding="utf-8"))
if len(sys.argv) > 1: sys.exit()
code = r'''
import json
from qgis.core import QgsProject, QgsGeometry, QgsFeature, QgsWkbTypes
p = json.load(open(PAY, encoding="utf-8")); l = QgsProject.instance().mapLayersByName("FMB_Fine_Tunned")[0]
assert not l.isEditable() and not l.isModified(), "the layer is being edited"; n0 = l.featureCount(); idx = {f.name(): i for i, f in enumerate(l.fields())}
feats = {}
for f in l.getFeatures():
    if "cut from the merged row" in str(f["topo_note"]): feats.setdefault((str(f["kide"]), round(float(f["area_sqm"]), 1)), []).append(f)
plan = []
for x in p["fix"]:
    g = QgsGeometry(); g.fromWkb(bytes.fromhex(x["wkb"])); c = [f for f in feats.get((x["kide"], round(x["area_old"], 1)), []) if f.geometry().isGeosEqual(g)]; assert len(c) == 1, (x["kide"], len(c)); plan.append((c[0].id(), x["attrs"]))
l.startEditing()
try:
    for fid, a in plan: assert l.changeAttributeValues(fid, {idx[c]: v for c, v in a.items()}), "attribute change refused"
    new = []
    for a in p["add"]:
        g = QgsGeometry(); g.fromWkb(bytes.fromhex(a["wkb"])); assert not g.isNull() and g.isGeosValid() and not g.isMultipart(), "a part is not a valid single polygon"
        f = QgsFeature(l.fields()); f.setGeometry(g)
        for c, v in a["attrs"].items(): f.setAttribute(idx[c], v)
        new.append(f)
    assert l.addFeatures(new), "parts refused"; ok = l.commitChanges(); err = [] if ok else l.commitErrors()[:5]
    if not ok: l.rollBack()
except Exception:
    l.rollBack(); raise
l.triggerRepaint(); print(json.dumps({"committed": ok, "errors": err, "features_before": n0, "features": l.featureCount()}))
'''.replace("PAY", repr(pf.replace("\\", "/")))
print(qb.json_result(code, timeout=900))
