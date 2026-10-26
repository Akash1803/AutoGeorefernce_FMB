"""Match every placed plot to its polygon on the sheet by POSITION, not by poly_id (ids of equal-sized plots swap between
conversions) and not by one rigid fit (hand-placed surveys are bent against their sheet): the sheet is warped onto the
layer with a thin-plate spline through the plots whose number is unique on both sides, then each plot takes the sheet
polygon it overlaps most. Labels = every machine read inside that sheet polygon. -> out/ledger4.csv, out/warp/<v>_<s>.gpkg"""
import warnings, sys, os, json; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.ops import unary_union
from scipy.interpolate import RBFInterpolator
from autogeoref import visible, fabric_check as fc
R = r"D:\Projects\Tambaram_Chengalpattu"
buf = unary_union(list(gpd.read_file(os.path.join(R, "Tambaram_Chengalpattu_Railway.gpkg"), layer="rail_buffer_30m").to_crs(32644).geometry))
clean = lambda x: "" if x is None or str(x).strip() in ("nan", "None", "NULL", "0", "") else str(x).strip()
AL = json.load(open("out/align2.json")); L2 = pd.read_csv("out/labels2.csv", dtype=str).fillna(""); L3 = pd.read_csv("out/labels3.csv", dtype=str).fillna("")
L2 = L2[(L2.ok == "True") & (L2.pid != "-1")]; num = L3.text.str.isdigit()
L3 = L3[(L3.pid != "-1") & (((~num) & (L3["isdim"] == "False")) | (num & (L3["strict"] == "False") & (L3.text.str.len() <= 3)))]
rows = []; os.makedirs("out/warp", exist_ok=True)
for V in ["35_04_074", "35_04_077"]:
    sh = visible.shown(V)
    for s in sorted(sh):
        S = gpd.read_file(f"out/fresh/{V}/{s}_parcels.geojson"); S["geometry"] = S.geometry.buffer(0)
        L = gpd.read_file(sh[s][0], layer=sh[s][1]); L = L.to_crs(32644) if L.crs.to_epsg() != 32644 else L; L = L[L.geometry.notna() & ~L.geometry.is_empty].copy()
        L["g"] = [max(fc._polys(g.buffer(0)), key=lambda q: q.area) for g in L.geometry]; L["lab"] = [clean(x) for x in L.plot_no]
        if len(L) < 2:
            r = L.iloc[0]; rows.append({"village_code": V, "survey": s, "fid": L.index[0], "layer_plot": "" if r.lab == s else r.lab, "single": True, "sheet_pid": int(S.poly_id.iloc[0]) if len(S) else -1, "w_iou": 1.0, "w_second": 0.0,
                                        "labels": "", "area": round(r.g.area, 1), "in_buffer": bool(r.g.intersects(buf)), "n_layer": 1, "n_sheet": len(S), "controls": 0, "map": "single"}); continue
        labs = {int(p): set() for p in S.poly_id}
        for p, l in zip(S.poly_id, S.plot_no):
            if clean(l): labs[int(p)].add(clean(l))
        for src in (L2, L3):
            x = src[(src.village_code == V) & (src.survey == s)]
            for p, t in zip(x.pid, x.text): labs[int(p)].add(t)
        a = AL[f"{V}|{s}"]; A0 = np.array([[a[0], a[1]], [a[2], a[3]]]); t0 = np.array([a[4], a[5]]); rig = lambda P: np.asarray(P) @ A0.T + t0
        # control pairs: a number that sits in exactly one sheet polygon and on exactly one placed plot
        cnt = {}
        for p, ls in labs.items():
            for l in ls: cnt.setdefault(l, []).append(p)
        lc = L.lab.value_counts(); P, Q = [], []
        for fid, r in L.iterrows():
            if r.lab and lc[r.lab] == 1 and len(cnt.get(r.lab, [])) == 1 and len(r.lab) >= 2:
                sp = S[S.poly_id == cnt[r.lab][0]].geometry.iloc[0]; P.append(sp.centroid.coords[0]); Q.append(r.g.centroid.coords[0])
        P, Q = np.array(P), np.array(Q); nctl = 0; maps = {"rigid": rig}
        if len(P) >= 4:
            res = np.hypot(*(rig(P) - Q).T); keep = res < max(6.0, 3 * np.median(res)); nctl = int(keep.sum())
            if nctl >= 4:
                X1 = np.c_[rig(P[keep]), np.ones(nctl)]; C, *_ = np.linalg.lstsq(X1, Q[keep], rcond=None)
                maps["affine"] = lambda X, C=C: np.c_[rig(X), np.ones(len(X))] @ C
            if nctl >= 12:
                tps = RBFInterpolator(rig(P[keep]), Q[keep] - rig(P[keep]), kernel="thin_plate_spline", smoothing=5.0)
                maps["warp"] = lambda X, tps=tps: rig(X) + tps(rig(X))
        best = None
        for nm, f in maps.items():
            SWk = [shapely.transform(g, lambda c: f(c)).buffer(0) for g in S.geometry]; tr = shapely.STRtree(SWk)
            sc = np.mean([max([r.g.intersection(SWk[k]).area / max(r.g.union(SWk[k]).area, 1e-9) for k in tr.query(r.g)] or [0]) for _, r in L.iterrows()])
            if best is None or sc > best[0] + 0.005: best = (sc, nm, SWk)
        how, SW = best[1], best[2]
        gpd.GeoDataFrame({"poly_id": S.poly_id.values, "labels": ["/".join(sorted(labs[int(p)])) for p in S.poly_id]}, geometry=SW, crs=32644).to_file(f"out/warp/{V[-3:]}_{s}.gpkg", driver="GPKG")
        tree = shapely.STRtree(SW)
        for fid, r in L.iterrows():
            ov = sorted([(r.g.intersection(SW[k]).area / max(r.g.union(SW[k]).area, 1e-9), k) for k in tree.query(r.g)], reverse=True)[:2]
            k = ov[0][1] if ov else -1; pid = int(S.poly_id.iloc[k]) if k >= 0 else -1
            rows.append({"village_code": V, "survey": s, "fid": fid, "layer_plot": r.lab, "single": False, "sheet_pid": pid, "w_iou": round(ov[0][0], 3) if ov else 0, "w_second": round(ov[1][0], 3) if len(ov) > 1 else 0,
                         "labels": "/".join(sorted(labs.get(pid, set()))), "area": round(r.g.area, 1), "in_buffer": bool(r.g.intersects(buf)), "n_layer": len(L), "n_sheet": len(S), "controls": nctl, "map": how})
X = pd.DataFrame(rows); X.index.name = "rid"; X.to_csv("out/ledger4.csv")
M = X[~X.single]; print("multi-plot rows", len(M), "| overlap with the matched sheet polygon: median %.3f, below 0.6: %d, below 0.4: %d" % (M.w_iou.median(), (M.w_iou < 0.6).sum(), (M.w_iou < 0.4).sum()))
d = M.groupby(["village_code", "survey", "sheet_pid"]).size(); print("sheet polygons claimed by more than one placed plot:", int((d > 1).sum()))
print(M[M.w_iou < 0.6][["village_code", "survey", "layer_plot", "labels", "w_iou", "w_second", "area", "in_buffer", "controls", "map"]].to_string()); print(M.drop_duplicates(["village_code", "survey"])["map"].value_counts().to_dict())
