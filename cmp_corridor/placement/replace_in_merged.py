"""Replace one survey in the merged layer by a georeferenced sheet (the parcel-by-parcel correction Akash chose):
- drop the CMP plots of that survey, insert the sheet's plots (red, source/status/pdf_check/placed_by/poly_id/merge_note filled),
- cut CMP parcels of other surveys where they overlap the new sheet (a CMP plot >= 95 % under it is removed and listed),
- fill gaps/slivers between the new sheet and CMP parcels into the CMP side (longest contact), never into a georeferenced sheet,
- report overlaps involving any georeferenced sheet (must be zero) and write the merged file (which must not be open in QGIS).
Usage: python replace_in_merged.py <merged.geojson> <village_code> <survey> <sheet.gpkg> "<placed_by note>" """
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
import pandas as pd
F, VC, SV, SHEET, NOTE = sys.argv[1:6]
assert ("CMP_my_Part" in F or "CMP_FMB_Fine_Tunned_Merged" in F) and os.path.exists(F) and os.path.exists(SHEET)
M = gpd.read_file(F).to_crs(32644)
# normalise: polygons only (a GeometryCollection left by make_valid has no boundary)
M["geometry"] = [as_poly(x) if (x is not None and x.geom_type not in ("Polygon", "MultiPolygon")) else x for x in M.geometry]
M["kide"] = M.kide.astype(str)
M["ours"] = M.source.astype(str).str.startswith("georef")
M["survey"] = M.kide.str.split("/").str[0]
M["village_code"] = M.village_code.astype(str)
log = []


def note(i, txt):
    M.loc[i, "merge_note"] = (str(M.loc[i, "merge_note"] or "") + "; " + txt).strip("; ")


old_cmp = M.index[(M.village_code == VC) & (M.survey == SV) & (~M.ours)]
old_geo = M.index[(M.village_code == VC) & (M.survey == SV) & (M.ours)]
old = old_cmp.append(old_geo)
log.append(f"{VC} {SV}: {len(old_cmp)} CMP plots and {len(old_geo)} earlier georeferenced plots removed from the merged layer")
src0 = gpd.read_file("C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2/merged_v1_backup.geojson")
src0["kide"] = src0.kide.astype(str); src0["village_code"] = src0.village_code.astype(str)
tpl = src0[(src0.village_code == VC) & (src0.kide.str.split("/").str[0] == SV)]
template = M.loc[old_cmp[0]] if len(old_cmp) else (tpl.iloc[0] if len(tpl) else None)
M = M.drop(old)
g = gpd.read_file(SHEET).to_crs(32644)
g["village_code"] = VC
g["survey"] = SV
g["ours"] = True
g["kide"] = g.kide.astype(str)
g["survey_no"] = SV
g["subdiv_no"] = g.subdiv_no.astype(str).replace({"-": ""})
g["unique_id"] = [f"{VC}_{k}" for k in g.kide]
for c in ("district", "taluk", "village"):
    g[c] = template[c] if template is not None else None
g["source"] = f"georeferenced sheet ({os.path.basename(SHEET)})"
g["status"] = "red - for review"
g["placed_by"] = NOTE
g["pdf_check"] = "numbers from the sketch vector; shape from the PDF line work"
g["merge_note"] = ""
# register columns by kide from the CMP rows just removed (same survey) when present
reg = [c for c in M.columns if c.startswith("tamilnilam_") or c.startswith("reginet_")]
lut = tpl.drop_duplicates("kide").set_index("kide") if len(tpl) else gpd.GeoDataFrame(columns=reg + ["kide"]).set_index("kide")
hits = 0
for c in reg:
    vals = []
    for k in g.kide:
        vals.append(lut.loc[k, c] if k in lut.index else None)
    g[c] = vals
hits = sum(1 for k in g.kide if k in lut.index)
log.append(f"register data carried over for {hits} of {len(g)} plots")
cols = [c for c in M.columns if c != "geometry"]
for c in cols:
    if c not in g.columns:
        g[c] = None
M = gpd.GeoDataFrame(pd.concat([M, g[cols + ["geometry"]]], ignore_index=True), geometry="geometry", crs=32644)
M["kide"] = M.kide.astype(str)
UN = make_valid(unary_union(list(g.geometry)))
# the new sheet yields to georeferenced sheets already in the layer (they were checked first)
others_geo = [x for x, o, s in zip(M.geometry, M.ours, M.survey) if o and s != SV and x is not None]
if others_geo:
    UG = make_valid(unary_union(others_geo))
    cut_new = []
    for i in M.index[M.ours & (M.survey == SV) & (M.village_code == VC)]:
        x = M.loc[i, "geometry"]
        a = x.intersection(UG).area
        if a > 0.5:
            y = as_poly(make_valid(x.difference(UG)))
            if y.is_empty or y.area < 0.5:
                cut_new.append((M.loc[i, "kide"], round(a), "entirely under a placed sheet - kept, flagged"))
                note(i, "lies entirely under a georeferenced sheet placed earlier - needs a look")
                continue
            M.loc[i, "geometry"] = y
            M.loc[i, "area_sqm"] = round(y.area, 1)
            M.loc[i, "area_acre"] = round(y.area / 4046.856, 4)
            note(i, f"cut by {a:.0f} m2 where it overlapped a georeferenced sheet placed earlier")
            cut_new.append((M.loc[i, "kide"], round(a), round(x.area), round(y.area)))
    log.append(f"{SV} cut where it met earlier georeferenced sheets: {cut_new}")
    UN = make_valid(unary_union([x for x, o, s in zip(M.geometry, M.ours, M.survey) if o and s == SV and x is not None]))
# cut other CMP parcels at the new sheet
removed, cuts = [], []
for i in list(M.index[(~M.ours)]):
    x = M.loc[i, "geometry"]
    if x is None or not x.intersects(UN):
        continue
    a = x.intersection(UN).area
    if a < 0.5:
        continue
    if a >= 0.95 * x.area:
        removed.append((M.loc[i, "kide"], round(x.area)))
        M = M.drop(i)
        continue
    y = as_poly(make_valid(x.difference(UN)))
    M.loc[i, "geometry"] = y
    M.loc[i, "area_sqm"] = round(y.area, 1)
    M.loc[i, "area_acre"] = round(y.area / 4046.856, 4)
    note(i, f"cut by {a:.0f} m2 where it overlapped georeferenced sheet {SV}")
    cuts.append((M.loc[i, "kide"], round(a), round(x.area), round(y.area)))
log.append(f"CMP parcels cut at {SV}: {cuts}")
log.append(f"CMP parcels removed (>= 95 % under {SV}): {removed}")
# any georeferenced sheet overlapping the new one? (should not happen; report)
UO_other = make_valid(unary_union([x for x, o, s in zip(M.geometry, M.ours, M.survey) if o and s != SV and x is not None]))
ovn = UN.intersection(UO_other).area
log.append(f"overlap of {SV} with other georeferenced sheets: {ovn:.1f} m2")
if os.environ.get('NOFILL') != '1':
    # gaps / slivers around the new sheet -> CMP side
    near = M[M.intersects(UN.buffer(60))]
    ALL = make_valid(unary_union([x for x in near.geometry if x is not None]))
    closed = ALL.buffer(3, join_style=2).buffer(-3, join_style=2)
    pieces = [p for p in getattr(make_valid(closed.difference(ALL)), "geoms", [closed.difference(ALL)]) if p.geom_type == "Polygon" and p.area > 1]
    for pp in getattr(ALL, "geoms", [ALL]):
        for r in pp.interiors:
            h = Polygon(r)
            if 1 < h.area < 500:
                pieces.append(h)
    filled, left = [], []
    for h in pieces:
        if h.distance(UN) > 0.1:
            continue
        cand = []
        for i in near.index:
            if i not in M.index or M.loc[i, "ours"]:
                continue
            gi = M.loc[i, "geometry"]
            bd = getattr(gi, "boundary", None) if gi is not None else None
            L = h.buffer(0.1).intersection(bd).length if bd is not None else 0
            if L > 0.5:
                cand.append((L, i))
        if not cand:
            left.append((round(h.area), round(h.centroid.x), round(h.centroid.y)))
            continue
        cand.sort(reverse=True)
        L, i = cand[0]
        gi = M.loc[i, "geometry"]
        y = as_poly(make_valid(unary_union([gi, h.buffer(0.02, join_style=2)])))
        y = as_poly(make_valid(y.difference(unary_union([z for j, z in zip(M.index, M.geometry) if j != i and z is not None and z.distance(h) < 1]))))
        M.loc[i, "geometry"] = y
        M.loc[i, "area_sqm"] = round(y.area, 1)
        M.loc[i, "area_acre"] = round(y.area / 4046.856, 4)
        note(i, f"stretched by {h.area:.0f} m2 to close the gap to georeferenced sheet {SV}")
        filled.append((M.loc[i, "kide"], round(h.area)))
    log.append(f"gaps filled into CMP parcels: {filled}")
    log.append(f"gaps left open (no CMP neighbour): {left}")
else:
    log.append('gaps left as they are (no gap filling, Akash 2026-10-23)')
M = M.reset_index(drop=True)
sidx = M.sindex
ov = []
for i, r in M.iterrows():
    if r.geometry is None:
        continue
    for j in sidx.query(r.geometry, predicate="intersects"):
        if j <= i:
            continue
        if not (r.ours or M.ours.iloc[j]):
            continue
        a = r.geometry.intersection(M.geometry.iloc[j]).area
        if a > 0.5:
            ov.append((r.kide, M.kide.iloc[j], round(a, 1)))
log.append(f"overlaps involving georeferenced sheets after the replacement: {ov}")
bad = [k for k, x in zip(M.kide, M.geometry) if x is None or x.is_empty or not x.is_valid]
log.append(f"invalid/empty: {bad}")
print("\n".join(log))
# final validity pass
M["geometry"] = [g if (g is None or g.is_valid) else as_poly(make_valid(g)) for g in M.geometry]
M.drop(columns=["ours", "survey"]).to_crs(4326).to_file(F, driver="GeoJSON")
print(f"written {F}: {len(M)} plots")
