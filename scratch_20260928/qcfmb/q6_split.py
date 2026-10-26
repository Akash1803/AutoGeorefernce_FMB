"""Dual QC, step 6: rows of the layer that merge several different plots of the sheet are cut back into those plots.
No outer boundary moves: the pieces are the sheet polygons as the pipeline placed them, clipped to the row, and together
they fill the row exactly. Each piece takes the number the second read found in its own sheet polygon, under the same
two-source rule as every other row. -> out/split.json, out/split_pieces.csv, out/visual2.csv"""
import os, re, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
from pyproj import Transformer
HERE = os.path.dirname(os.path.abspath(__file__)); R = r"D:\Projects\Tambaram_Chengalpattu"
G4 = gpd.read_file(os.path.join(HERE, "out", "layer_before.geojson")); G = G4.to_crs(32644)
D = pd.read_csv(os.path.join(HERE, "out", "plan.csv"), dtype=str).fillna(""); K = pd.read_pickle(os.path.join(HERE, "out", "link_geom.pkl")); K = K[K.wkb.astype(str).str.len() > 0]; K["share"] = K.share.astype(float); K["row"] = K.row.astype(int)
L = pd.read_csv(os.path.join(HERE, "out", "labels5_all.csv"), dtype=str).fillna("")
good = L[(L.plotno == "True") & (L.pid != "-1")].groupby(["village_code", "sheet", "pid"]).text.apply(lambda s: sorted(set(s))).to_dict()
T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3); T["sqm"] = T.extent_hectares.astype(float) * 10000 + T.extent_ares.astype(float) * 100; T = T[T.fmb_survey_no == T.clip_survey_no]
reg = {k: dict(zip(g.sub_division, g.sqm)) for k, g in T.groupby(["vc", "fmb_survey_no"])}
to4326 = Transformer.from_crs(32644, 4326, always_xy=True)
polys = lambda g: [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if p.geom_type == "Polygon" and p.area > 0]
out, rows, bad = [], [], []
for r in D[D.merged == "True"].itertuples():
    i = int(r.row); GR = G.geometry[i].buffer(0); k = K[(K.row == i) & (K.share >= 0.5)]
    P = [(str(p), shapely.snap(shapely.from_wkb(bytes.fromhex(w)).buffer(0), GR, 0.05).buffer(0)) for p, w in zip(k.pid, k.wkb)]
    P.sort(key=lambda t: -t[1].area); taken = None; Q = []
    for pid, g in P:
        q = g.intersection(GR); q = q if taken is None else q.difference(taken); q = unary_union(polys(q)) if polys(q) else None
        if q is None or q.area < 0.01: continue
        Q.append([pid, q]); taken = q if taken is None else unary_union([taken, q])
    if len(Q) < 2: bad.append((i, r.kide, "fewer than two pieces")); continue
    # what the placed polygons leave of the row goes to the piece it borders most; loose bits of a piece go to their neighbour
    for _pass in range(3):
        rest = polys(GR.difference(unary_union([q for _p, q in Q])))
        for j, (pid, q) in enumerate(Q):
            pp = polys(q)
            if len(pp) > 1: pp.sort(key=lambda x: -x.area); Q[j][1] = pp[0]; rest += pp[1:]
        if not rest: break
        for b in rest:
            sh = [(b.buffer(0.002).intersection(q.boundary).length, j) for j, (_p, q) in enumerate(Q)]; l, j = max(sh)
            if l <= 0: l, j = max((-(b.distance(q)), j) for j, (_p, q) in enumerate(Q))
            Q[j][1] = unary_union([Q[j][1], b.buffer(0)]).buffer(0)
    tot = sum(q.area for _p, q in Q); ov = sum(Q[a][1].intersection(Q[b][1]).area for a in range(len(Q)) for b in range(a + 1, len(Q)))
    if abs(tot - GR.area) > 0.05 or ov > 0.01 or any(not q.is_valid for _p, q in Q): bad.append((i, r.kide, "pieces do not fill the row: %.2f vs %.2f, overlap %.3f" % (tot, GR.area, ov))); continue
    RG = reg.get((r.village_code, r.sheet)); pieces = []
    for pid, q in Q:
        c = [t for t in good.get((r.village_code, r.sheet, pid), []) if t != r.sheet and t != re.sub(r"[A-Z]+$", "", r.sheet)]
        inreg = lambda t: RG is not None and t in RG; fits = lambda t: inreg(t) and RG[t] > 0 and (0.5 <= q.area / RG[t] <= 2.0 or abs(q.area - RG[t]) <= 50)
        strong = [t for t in c if len(t) >= 2 or inreg(t)]; f = [t for t in strong if fits(t)]
        if len(f) == 1: n, s, v = f[0], "second read + A-register", False
        elif len(strong) >= 1: n, s, v = "", "to be looked at on the PDF (second read: %s)" % "/".join(strong), True
        else: n, s, v = "", "no number printed on the sheet (both reads)", False
        pieces.append({"pid": pid, "sub": n, "source": s, "visual": v, "geom": q, "second_read": "/".join(c)})
    # true pieces of one subdivision that touch become one plot again
    merged = []
    for p in pieces:
        m = next((x for x in merged if p["sub"] and x["sub"] == p["sub"] and x["geom"].buffer(0.01).intersects(p["geom"])), None)
        if m: m["geom"] = unary_union([m["geom"], p["geom"]]).buffer(0); m["pid"] += "/" + p["pid"]
        else: merged.append(dict(p))
    # a number set on two pieces that do not touch is looked at
    cnt = pd.Series([p["sub"] for p in merged if p["sub"]]).value_counts()
    for p in merged:
        if p["sub"] and cnt[p["sub"]] > 1: p["visual"] = True; p["source"] = "to be looked at on the PDF (number %s would repeat)" % p["sub"]; p["second_read"] = p["second_read"] or p["sub"]; p["sub"] = ""
    # the same cut on the layer's own coordinates (lon/lat, 8 decimals): the outer line of the row is kept vertex for vertex,
    # only the lines between the pieces are new
    row4 = G4.geometry[i]; t4 = lambda g: shapely.make_valid(shapely.transform(g, lambda c: np.column_stack(to4326.transform(c[:, 0], c[:, 1]))))
    P4 = [shapely.snap(t4(p["geom"]), row4, 3e-7).buffer(0) for p in merged]; order = sorted(range(len(P4)), key=lambda k: -P4[k].area); F4 = [None] * len(P4); taken = None
    for k in order:
        q = P4[k].intersection(row4); q = q if taken is None else q.difference(taken); pp = [x for x in shapely.get_parts(q) if x.geom_type == "Polygon" and x.area > 1e-14]
        F4[k] = shapely.unary_union(pp) if pp else None; taken = F4[k] if taken is None else (taken if F4[k] is None else shapely.unary_union([taken, F4[k]]))
    if any(f is None for f in F4): bad.append((i, r.kide, "a piece vanished on the layer's coordinates")); rows = [x for x in rows if x["row"] != i]; continue
    for _pass in range(3):
        rest = [x for x in shapely.get_parts(row4.difference(shapely.unary_union(F4))) if x.geom_type == "Polygon" and x.area > 0]
        if not rest: break
        for b_ in rest:
            sh = [(b_.buffer(2e-8).intersection(f.boundary).length, k) for k, f in enumerate(F4)]; l_, k = max(sh)
            if l_ <= 0: k = int(np.argmin([b_.distance(f) for f in F4]))
            F4[k] = shapely.unary_union([F4[k], b_])
    F4 = [shapely.set_precision(f, 1e-8) for f in F4]; U4 = shapely.unary_union(F4)
    ovl = sum(F4[a].intersection(F4[b]).area for a in range(len(F4)) for b in range(a + 1, len(F4)))
    if any(f.is_empty or not f.is_valid for f in F4) or U4.symmetric_difference(row4).area > 2e-14 or ovl > 2e-14:
        bad.append((i, r.kide, "the pieces do not tile the row: difference %.2e, overlap %.2e (square degrees)" % (U4.symmetric_difference(row4).area, ovl))); rows = [x for x in rows if x["row"] != i]; continue
    # one feature, one polygon: a piece left in several parts keeps its largest part; a part under 5 m2 joins the piece it lies
    # along, a larger part stands as its own plot under the same number
    sqm = lambda g: float(gpd.GeoSeries([g], crs=4326).to_crs(32644).area.iloc[0]); extra = []
    for k in range(len(F4)):
        parts = sorted(shapely.get_parts(F4[k]), key=lambda x: -x.area)
        if len(parts) == 1: continue
        F4[k] = parts[0]
        for b_ in parts[1:]:
            if sqm(b_) < 5.0:
                sh = [(b_.buffer(2e-8).intersection(F4[j].boundary).length if j != k else -1, j) for j in range(len(F4))]; l_, j = max(sh); u_ = shapely.unary_union([F4[j], b_]) if l_ > 0 else None
                if u_ is not None and u_.geom_type == "Polygon": F4[j] = shapely.set_precision(u_, 1e-8); continue
            extra.append((k, b_))
    for k, b_ in extra: merged.append(dict(merged[k])); F4.append(b_)
    U4 = shapely.unary_union(F4); ovl = sum(F4[a].intersection(F4[b]).area for a in range(len(F4)) for b in range(a + 1, len(F4)))
    if any(f.geom_type != "Polygon" or not f.is_valid for f in F4) or U4.symmetric_difference(row4).area > 2e-14 or ovl > 2e-14:
        bad.append((i, r.kide, "the pieces do not tile the row after the parts were settled")); rows = [x for x in rows if x["row"] != i]; continue
    for p, g4 in zip(merged, F4):
        am = sqm(g4); p["area4"] = am
        rows.append({"row": i, "village_code": r.village_code, "village": r.village, "sheet": r.sheet, "old_kide": r.kide, "pid": p["pid"], "sub": p["sub"], "source": p["source"], "visual": p["visual"], "second_read": p["second_read"], "area": round(p["area4"], 1), "area_utm": round(p["geom"].area, 1), "wkb": g4.wkb_hex})
    out.append(i)
S = pd.DataFrame(rows); S.to_csv(os.path.join(HERE, "out", "split_pieces.csv"), index=False)
print("merged rows", int((D.merged == "True").sum()), "| split", len(out), "| not split", len(bad)); [print("   ", b) for b in bad]
g_ = [shapely.from_wkb(bytes.fromhex(w)) for w in S.wkb]; print("pieces valid on the layer's grid:", sum(x.is_valid for x in g_), "of", len(g_), "| multi-part:", sum(x.geom_type != "Polygon" for x in g_), "| largest area change against the metric piece: %.2f m2" % (S.area.astype(float) - S.area_utm.astype(float)).abs().max())
print("pieces", len(S), "| numbered by two sources", int((S["sub"] != "").sum()), "| to look at on the PDF", int(S.visual.sum()), "| no number printed", int(((S["sub"] == "") & ~S.visual).sum()))
V = S[S.visual].copy(); V["pids"] = V.pid; V["row"] = ["%d_%s" % (a, b.split("/")[0]) for a, b in zip(V.row, V.pid)]; V = V.drop_duplicates("row"); V.to_csv(os.path.join(HERE, "out", "visual2.csv"), index=False); S["piece"] = ["%d_%s" % (a, b.split("/")[0]) for a, b in zip(S.row, S.pid)]; S.to_csv(os.path.join(HERE, "out", "split_pieces.csv"), index=False)
