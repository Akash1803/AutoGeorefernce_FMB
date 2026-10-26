"""Place the 124 subdivisions that are drawn on the PDF but missing from CMP_two_Parts into the unnumbered '124' lumps:
for each lump, a LOCAL rigid pose of the sheet from the numbered plots within 150 m (label pairs, outliers dropped); every missing
sheet plot lying mostly in the lump is placed with that pose and cut to the lump (existing plots never move); the rest of the lump
stays '124' (road). Writes fix_two_fill124.json (lump geometry change + added plots) + fill124.png"""
import sys, os, math, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd, geopandas as gpd
from shapely import make_valid
from shapely.ops import unary_union
from shapely.affinity import affine_transform
exec(open("georef124.py", encoding="utf-8").read().split("L = gpd.read_file")[0])          # sheet plots S (sheet frame) + labels
norm = lambda s: "".join(ch for ch in str(s).upper() if ch.isalnum())
L = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); L["kide"] = L.kide.astype(str)
cur = L[(L.village_code.astype(str) == VC) & (L.kide.str.split("/").str[0] == SV)].copy(); cur["s"] = cur.subdiv_no.astype(str).map(norm).replace({"NONE": "", "NAN": ""})
def kabsch(P, Q):
    cp, cq = P.mean(0), Q.mean(0); Hm = (P - cp).T @ (Q - cq); U, Sg, Vt = np.linalg.svd(Hm); Rm = Vt.T @ U.T
    if np.linalg.det(Rm) < 0: Vt[-1] *= -1; Rm = Vt.T @ U.T
    return Rm, cq - Rm @ cp
pairs = []
for t in set(S.subdiv_no) & set(cur.s) - {""}:
    a, b = S[S.subdiv_no == t], cur[cur.s == t]
    if len(a) == 1 and len(b) == 1: pairs.append((t, np.array(a.geometry.iloc[0].representative_point().coords[0]), np.array(b.geometry.iloc[0].representative_point().coords[0])))
missing = S[(S.subdiv_no != "") & ~S.subdiv_no.isin(set(cur.s))]
lumps = cur[cur.kide == SV]
out = {"change": [], "add": [], "delete": []}; rep = []; placed_all = []
tmpl = cur[cur.kide != SV].iloc[0]
for li, lump in lumps.iterrows():
    LG = make_valid(lump.geometry); c = np.array(LG.centroid.coords[0])
    near = [p for p in pairs if np.hypot(*(p[2] - c)) < 150]
    if len(near) < 4: near = sorted(pairs, key=lambda p: np.hypot(*(p[2] - c)))[:8]
    P0 = np.array([p[1] for p in near]); Q0 = np.array([p[2] for p in near]); keep = np.ones(len(near), bool)
    for _ in range(5):
        Rm, t = kabsch(P0[keep], Q0[keep]); res = np.hypot(*((P0 @ Rm.T + t) - Q0).T); keep = res < max(2.0, 2.5 * np.median(res[keep]))
    Rm, t = kabsch(P0[keep], Q0[keep]); res = np.hypot(*((P0 @ Rm.T + t) - Q0).T)
    Mx = [Rm[0, 0], Rm[0, 1], Rm[1, 0], Rm[1, 1], t[0], t[1]]
    new = []
    for r in missing.itertuples():
        g = affine_transform(r.geometry, Mx)
        if g.intersection(LG).area < 0.6 * g.area: continue
        q = make_valid(g.intersection(LG)); q = max([z for z in getattr(q, "geoms", [q]) if z.geom_type == "Polygon"], key=lambda z: z.area)
        if q.area < 5: continue
        new.append((r.subdiv_no, q, g.area))
    if not new: continue
    rest = make_valid(LG.difference(unary_union([q.buffer(0.001, join_style=2) for _, q, _ in new])))
    parts = sorted([z for z in getattr(rest, "geoms", [rest]) if z.geom_type == "Polygon"], key=lambda z: -z.area)
    rep.append({"lump_m2": round(LG.area), "pose_pairs": int(keep.sum()), "pose_residual_m": round(float(np.median(res[keep])), 2), "placed": [(s, round(q.area), round(a)) for s, q, a in new],
                "road_left_m2": round(sum(z.area for z in parts)), "road_pieces": len(parts)})
    # the lump keeps the road remainder (largest piece); other remainder pieces stay unnumbered 124 too (added as their own rows)
    out["change"].append({"vc": VC, "old_kide": SV, "old_pt": list(LG.representative_point().coords[0]), "wkt": parts[0].wkt if parts else LG.wkt,
                          "attrs": {"area_sqm": round(parts[0].area, 1) if parts else round(LG.area, 1), "area_acre": round((parts[0].area if parts else LG.area) / 4046.856, 4)}})
    for z in parts[1:]:
        if z.area < 1: continue
        a = {k: (None if str(v) in ("nan", "None") else v) for k, v in lump.drop(labels=["geometry", "s"]).to_dict().items()}; a.update({"area_sqm": round(z.area, 1), "area_acre": round(z.area / 4046.856, 4)})
        out["add"].append({"wkt": z.wkt, "attrs": a})
    for s, q, a in new:
        at = {k: (None if str(v) in ("nan", "None") else v) for k, v in tmpl.drop(labels=["geometry", "s"]).to_dict().items()}
        for k in list(at):
            if k.startswith("tamilnilam_") or k.startswith("reginet_"): at[k] = None
        at.update({"kide": f"{SV}/{s}", "survey_no": SV, "subdiv_no": s, "unique_id": f"{VC}_{SV}/{s}", "area_sqm": round(q.area, 1), "area_acre": round(q.area / 4046.856, 4),
                   "source": "georeferenced from the 124 PDF (plot missing in the vector conversion)", "status": "red - for review", "review": "rigid placement - for review"})
        out["add"].append({"wkt": q.wkt, "attrs": at}); placed_all.append((s, q))
for r in rep: print(r)
json.dump(out, open("fix_two_fill124.json", "w"), default=str); print("changes", len(out["change"]), "adds", len(out["add"]))
gpd.GeoDataFrame({"s": [s for s, _ in placed_all]}, geometry=[q for _, q in placed_all], crs=32644).to_file("fill124_preview.gpkg", driver="GPKG")
