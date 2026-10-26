"""Dual QC of FMB_Fine_Tunned, step 1: second, independent read of every sheet (read3: one character at a time, A-Z)
laid over the polygons the layer was built from (FMB_Vector/<village>/<sheet>_parcels.geojson).
-> out/labels_all.csv (every label read), out/polys_all.csv (one row per sheet polygon: first read, second read)"""
import sys, os, glob, warnings, traceback; warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import pandas as pd, geopandas as gpd, shapely
from multiprocessing import Pool
R = r"D:\Projects\Tambaram_Chengalpattu"
clean = lambda x: "" if x is None or str(x).strip().lower() in ("nan", "none", "null", "0", "") else str(x).strip()

def one(job):
    v, s, scale = job
    try:
        import read3
        S = gpd.read_file(os.path.join(R, "FMB_Vector", v, s + "_parcels.geojson")); S["geometry"] = S.geometry.buffer(0); S["l1"] = [clean(x) for x in S.plot_no]
        pdf = os.path.join(R, "FMB_Sketches", v, s + ".pdf")
        if not os.path.exists(pdf): return v, s, "no pdf", [], []
        labs, hL = read3.read(pdf, scale, known=[l for l in S.l1 if l])
        L = []
        for a in labs:
            p = shapely.Point(a["x"], a["y"]); hit = [k for k, q in enumerate(S.geometry) if q.contains(p)]
            k = min(hit, key=lambda k: S.geometry.iloc[k].area) if hit else -1
            L.append({"village_code": v, "sheet": s, "text": a["text"], "worst": a["worst"], "isdim": a["isdim"], "strict": a["strict"], "x": round(a["x"], 2), "y": round(a["y"], 2), "pid": int(S.poly_id.iloc[k]) if k >= 0 else -1, "cap_pt": round(hL, 2) if hL else ""})
        P = [{"village_code": v, "sheet": s, "pid": int(r.poly_id), "l1": r.l1, "alt1": clean(getattr(r, "alt_labels", "")), "area": round(r.geometry.area, 1), "n_sheet": len(S),
              "cx": round(r.geometry.representative_point().x, 2), "cy": round(r.geometry.representative_point().y, 2)} for r in S.itertuples()]
        return v, s, "ok" if hL else "no label font found", L, P
    except Exception as e:
        return v, s, "error: " + repr(e)[:200], [], []

if __name__ == "__main__":
    scales = {(r.village_code, r.survey_no): int(r.scale) for r in pd.read_csv(os.path.join(R, "FMB_Sketches", "sheet_scales.csv"), dtype=str).itertuples()}
    jobs = []
    for f in sorted(glob.glob(os.path.join(R, "FMB_Vector", "*", "*_parcels.geojson"))):
        v = os.path.basename(os.path.dirname(f)); s = os.path.basename(f)[:-len("_parcels.geojson")]; jobs.append((v, s, scales.get((v, s), 500)))
    print("sheets", len(jobs), flush=True); LL, PP, ST = [], [], []
    with Pool(6) as pool:
        for k, (v, s, st, L, P) in enumerate(pool.imap_unordered(one, jobs, chunksize=2)):
            LL += L; PP += P; ST.append({"village_code": v, "sheet": s, "state": st, "labels": len(L), "polys": len(P)})
            if k % 50 == 0: print(k, v, s, st, flush=True)
    pd.DataFrame(LL).to_csv(os.path.join(HERE, "out", "labels_all.csv"), index=False); pd.DataFrame(PP).to_csv(os.path.join(HERE, "out", "polys_all.csv"), index=False); S = pd.DataFrame(ST); S.to_csv(os.path.join(HERE, "out", "sheets_all.csv"), index=False)
    print("labels", len(LL), "| polygons", len(PP), "|", S.state.value_counts().to_dict())
