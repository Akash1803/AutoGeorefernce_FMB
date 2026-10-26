"""Dual QC of FMB_Fine_Tunned, step 2: which sheet polygon is which layer row.
Every sheet is put on its own rows of the layer the way the pipeline did it (rigid pose on the rows' outline, then the
conform warp), so each sheet polygon lands on the row made from it; the link is the row it overlaps most.
Checked with the rows whose number is unique on both sides. -> out/link_all.csv (one row per sheet polygon)"""
import sys, os, glob, warnings, logging; warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely import affinity
from shapely.geometry import Point
from shapely.ops import unary_union
from multiprocessing import Pool
R = r"D:\Projects\Tambaram_Chengalpattu"; LAYER = os.path.join(HERE, "out", "layer_before.geojson") if len(sys.argv) < 2 else sys.argv[1]
clean = lambda x: "" if x is None or str(x).strip().lower() in ("nan", "none", "null", "0", "") else str(x).strip()

def one(job):
    v, unit, rows = job              # rows: [(layer index, subdiv, wkb)]
    from autogeoref import fmb_on_base as fb, fmb_finish as ff
    out = []
    try:
        f = os.path.join(R, "FMB_Vector", v, unit + "_parcels.geojson")
        if not os.path.exists(f): return [{"village_code": v, "sheet": unit, "pid": -1, "row": i, "state": "no sheet for this unit"} for i, _s, _w in rows]
        S = gpd.read_file(f); S["geometry"] = [fb._valid(g) for g in S.geometry]; S = S[S.geometry.notna()]; S["l1"] = [clean(x) for x in S.plot_no]
        L = [(i, s, fb._valid(shapely.from_wkb(w))) for i, s, w in rows]; target = fb._valid(unary_union([g for _i, _s, g in L])); outline = unary_union(list(S.geometry))
        pose = fb.place(outline, target)
        if pose is None: return [{"village_code": v, "sheet": unit, "pid": int(p), "row": -1, "state": "no pose"} for p in S.poly_id]
        best = None
        for turn in (0.0, 180.0):
            P = dict(pose); P["angle"] = pose["angle"] + turn
            placed = [{"geometry": fb.apply_pose(g, P, target), "pid": int(p), "l1": l} for g, p, l in zip(S.geometry, S.poly_id, S.l1)]
            placed = [r for r in placed if r["geometry"] is not None and not r["geometry"].is_empty]
            try: ff.conform_survey(placed, target)
            except Exception: pass
            placed = [r for r in placed if r["geometry"] is not None and not r["geometry"].is_empty]
            tr = shapely.STRtree([g for _i, _s, g in L]); tot = 0.0; res = []
            for r in placed:
                c = sorted([(r["geometry"].intersection(L[k][2]).area, k) for k in tr.query(r["geometry"])], reverse=True)[:2]
                a = r["geometry"].area or 1e-9; k = c[0][1] if c and c[0][0] > 0 else -1
                res.append((r, k, c[0][0] / a if c else 0.0, c[1][0] / a if len(c) > 1 else 0.0)); tot += c[0][0] if c else 0.0
            # agreement of the numbers that are unique on both sides
            sub = pd.Series([s for _i, s, _g in L]); uq = set(sub[sub != ""].value_counts()[lambda x: x == 1].index) & set(S.l1[S.l1 != ""].value_counts()[lambda x: x == 1].index)
            agree = sum(1 for r, k, a1, a2 in res if r["l1"] in uq and k >= 0 and L[k][1] == r["l1"]); score = (agree, tot)
            if best is None or score > best[0]: best = (score, res, len(uq), agree, turn)
        # where the pose does not put the numbers on their rows: carry the sheet over through the plots whose number is unique on both sides
        if best[2] >= 3 and best[3] < 0.9 * best[2] or (best[2] >= 3 and sum(1 for r_, k_, a1_, a2_ in best[1] if a1_ < 0.5) > 0.15 * len(best[1])):
            from scipy.interpolate import RBFInterpolator
            sub = pd.Series([s_ for _i, s_, _g in L]); uq = set(sub[sub != ""].value_counts()[lambda x: x == 1].index) & set(S.l1[S.l1 != ""].value_counts()[lambda x: x == 1].index)
            P_ = np.array([S.geometry[S.l1 == u].iloc[0].representative_point().coords[0] for u in sorted(uq)]); Q_ = np.array([[g for _i, s_, g in L if s_ == u][0].representative_point().coords[0] for u in sorted(uq)])
            maps = {}
            X1 = np.c_[P_, np.ones(len(P_))]; C_, *_ = np.linalg.lstsq(X1, Q_, rcond=None); res0 = np.hypot(*((X1 @ C_) - Q_).T); keep = res0 <= max(8.0, 3 * np.median(res0))
            if keep.sum() >= 3:
                C2, *_ = np.linalg.lstsq(X1[keep], Q_[keep], rcond=None); maps["affine"] = lambda X, C2=C2: np.c_[X, np.ones(len(X))] @ C2
                if keep.sum() >= 12:
                    aff = maps["affine"]; tps = RBFInterpolator(aff(P_[keep]), Q_[keep] - aff(P_[keep]), kernel="thin_plate_spline", smoothing=5.0); maps["warp"] = lambda X, tps=tps, aff=aff: aff(X) + tps(aff(X))
            for nm, f in maps.items():
                tr = shapely.STRtree([g for _i, _s, g in L]); tot = 0.0; res = []
                for g, p, l in zip(S.geometry, S.poly_id, S.l1):
                    q = fb._valid(shapely.transform(g, lambda c: f(c)))
                    if q is None or q.is_empty: continue
                    c = sorted([(q.intersection(L[k][2]).area, k) for k in tr.query(q)], reverse=True)[:2]; a = q.area or 1e-9; k = c[0][1] if c and c[0][0] > 0 else -1
                    res.append(({"geometry": q, "pid": int(p), "l1": l}, k, min(c[0][0] / a, 1.0) if c else 0.0, c[1][0] / a if len(c) > 1 else 0.0)); tot += c[0][0] if c else 0.0
                agree = sum(1 for r_, k, a1, a2 in res if r_["l1"] in uq and k >= 0 and L[k][1] == r_["l1"])
                if (agree, tot) > best[0]: best = ((agree, tot), res, len(uq), agree, nm)
        _sc, res, nuq, agree, turn = best
        for r, k, a1, a2 in res:
            out.append({"village_code": v, "sheet": unit, "pid": r["pid"], "l1": r["l1"], "row": L[k][0] if k >= 0 else -1, "row_subdiv": L[k][1] if k >= 0 else "", "share": round(a1, 3), "second": round(a2, 3),
                        "wkb": r["geometry"].wkb_hex if os.environ.get("DUMP") else "", "placed_area": round(r["geometry"].area, 1), "row_area": round(L[k][2].area, 1) if k >= 0 else "", "controls": nuq, "controls_ok": agree, "turn": turn, "state": "ok"})
    except Exception as e:
        out.append({"village_code": v, "sheet": unit, "pid": -1, "row": -1, "state": "error: " + repr(e)[:200]})
    return out

if __name__ == "__main__":
    G = gpd.read_file(LAYER).to_crs(32644); G["unit"] = G.kide.astype(str).str.split("/").str[0]; G["sub"] = [clean(x) for x in G.subdiv_no]
    # the sheet a row was made from: the unit its plot_no carries (10A in "10A/2", or "10A" alone for an unnumbered plot of sheet 10A), else the front of kide
    import re
    def home(v, p, k, sv):
        p = clean(p); k = str(k).strip(); sv = str(sv).strip()
        if "/" in p: return p.split("/")[0].strip()
        if p and re.fullmatch(re.escape(sv) + r"[A-Z]{1,2}", p) and os.path.exists(os.path.join(R, "FMB_Vector", v, p + "_parcels.geojson")): return p       # an unnumbered plot of unit sheet 10A
        return k.split("/")[0].strip()
    G["home"] = [home(v, p, k, sv) for v, p, k, sv in zip(G.village_code, G.plot_no, G.kide, G.survey_no)]
    G["sub"] = ["" if s_ == h else s_ for s_, h in zip(G["sub"], G.home)]
    print("rows whose sheet is not the front of their kide:", int((G.home != G.unit).sum()), flush=True)
    G[["village_code", "kide", "unit", "home", "sub"]].to_csv(os.path.join(HERE, "out", "row_home.csv"))
    jobs = [(v, u, [(int(i), s, g.wkb) for i, s, g in zip(g_.index, g_["sub"], g_.geometry)]) for (v, u), g_ in G.groupby(["village_code", "home"])]
    print("layer rows", len(G), "| units", len(jobs), flush=True); rows = []
    with Pool(6) as pool:
        for k, r in enumerate(pool.imap_unordered(one, jobs, chunksize=2)):
            rows += r
            if k % 60 == 0: print(k, r[0]["village_code"], r[0]["sheet"], r[0]["state"], flush=True)
    X = pd.DataFrame(rows)
    if os.environ.get("DUMP"): X[["village_code", "sheet", "pid", "row", "share", "wkb"]].to_pickle(os.path.join(HERE, "out", "link_geom.pkl")); X = X.drop(columns="wkb")
    X.to_csv(os.path.join(HERE, "out", "link_all.csv"), index=False)
    ok = X[X.state == "ok"]; print("sheet polygons", len(X), "|", X.state.str.split(":").str[0].value_counts().to_dict())
    print("share of the polygon on its row: median %.3f | below 0.8: %d | below 0.5: %d" % (ok.share.median(), (ok.share < 0.8).sum(), (ok.share < 0.5).sum()))
    c = ok.drop_duplicates(["village_code", "sheet"]); print("numbers unique on both sides:", int(c.controls.sum()), "| landed on the row of the same number:", int(c.controls_ok.sum()))
    lab = ok[(ok.l1 != "")]; print("numbered polygons on a row with the same number: %d of %d" % ((lab.l1 == lab.row_subdiv).sum(), len(lab)))
    linked = set(ok.row[ok.row >= 0]); print("layer rows that received a sheet polygon:", len(linked), "of", len(G), "| rows without:", len(G) - len(linked))
