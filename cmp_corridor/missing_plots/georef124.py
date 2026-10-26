"""Kattour 124: plots of the PDF sheet (all solid black line work polygonized) + their printed numbers (read5, register-checked),
posed rigidly onto the numbered 124 plots already in CMP_two_Parts (label pairs, least squares + outline refine).
Reports which register subdivisions exist on the sheet but not in the layer, and where they land (inside the unnumbered 124 lumps?).
-> k124_sheet.gpkg (posed sheet plots with labels), printed summary"""
import sys, os, math, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "qc", "r5")); sys.path.insert(1, r"D:\code\FMB_to_GeoJSON")
import numpy as np, pandas as pd, geopandas as gpd, fitz, read5
from shapely import make_valid
from shapely.geometry import LineString, Point
from shapely.ops import unary_union, polygonize
from shapely.affinity import affine_transform
PT = 72 / 25.4 * 2
VC, SV = "35_11_117", "124"
PDF = f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/{VC}/{SV}.pdf"
norm = lambda s: "".join(ch for ch in str(s).upper() if ch.isalnum())
pg = fitz.open(PDF)[0]; W, H = pg.rect.width, pg.rect.height
lines = []
for d in pg.get_drawings():
    w = round(d.get("width") or 0, 2); col = tuple(d.get("color") or (0, 0, 0))
    if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1": continue
    for it in d["items"]:
        if it[0] == "l":
            x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
            if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
polys = [p for p in polygonize(unary_union(lines)) if p.area > 0.5]
T = pd.read_excel("D:/Projects/CMP_Roads_LPS/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
REG = {norm(x) for x in T[(T.vc == VC) & (T.fmb_survey_no == SV)].sub_division} - {""}
labs, info = read5.read(PDF, 500); ok = {id(a) for a in read5.plot_numbers(labs)}
nums = [(norm(a["text"]), Point(a["px"] / PT, -a["py"] / PT)) for a in labs if id(a) in ok]
nums_reg = [(t, p) for t, p in nums if t in REG]
print(f"sheet plots {len(polys)} | numbers read {len(nums)} | of them in the register {len(nums_reg)} | register size {len(REG)}")
# label each sheet plot: number inside (the one nearest its centre), else number within 3 m beside it (thin plots)
lab = []
for q in polys:
    inside = [(p.distance(q.representative_point()), t) for t, p in nums_reg if q.contains(p)]
    if inside: lab.append(min(inside)[1]); continue
    near = [(q.distance(p), t) for t, p in nums_reg if q.distance(p) <= 3.0]
    lab.append(min(near)[1] if near else "")
S = gpd.GeoDataFrame({"subdiv_no": lab}, geometry=polys)
dupl = S[S.subdiv_no != ""].subdiv_no.value_counts(); print("sheet plots with a number:", int((S.subdiv_no != "").sum()), "| numbers used twice:", dupl[dupl > 1].to_dict())
L = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); L["kide"] = L.kide.astype(str)
cur = L[(L.village_code.astype(str) == VC) & (L.kide.str.split("/").str[0] == SV)].copy(); cur["s"] = cur.subdiv_no.astype(str).map(norm).replace({"NONE": "", "NAN": ""})
# label pairs (unique on both sides) -> rigid fit (rotation + translation, true scale)
pairs = []
for t in set(S.subdiv_no) & set(cur.s) - {""}:
    a, b = S[S.subdiv_no == t], cur[cur.s == t]
    if len(a) == 1 and len(b) == 1: pairs.append((t, a.geometry.iloc[0].representative_point().coords[0], b.geometry.iloc[0].representative_point().coords[0], a.geometry.iloc[0].area, b.geometry.iloc[0].area))
P0 = np.array([p[1] for p in pairs]); Q0 = np.array([p[2] for p in pairs])
def kabsch(P, Q):
    cp, cq = P.mean(0), Q.mean(0); Hm = (P - cp).T @ (Q - cq); U, Sg, Vt = np.linalg.svd(Hm); Rm = Vt.T @ U.T
    if np.linalg.det(Rm) < 0: Vt[-1] *= -1; Rm = Vt.T @ U.T
    return Rm, cq - Rm @ cp
keep = np.ones(len(pairs), bool)
for it in range(5):                                   # drop pairs that disagree (mislabelled plots) and refit
    Rm, t = kabsch(P0[keep], Q0[keep]); res = np.hypot(*((P0 @ Rm.T + t) - Q0).T)
    lim = max(3.0, 3 * np.median(res[keep])); keep = res < lim
Rm, t = kabsch(P0[keep], Q0[keep]); res = np.hypot(*((P0 @ Rm.T + t) - Q0).T)
Mx = [Rm[0, 0], Rm[0, 1], Rm[1, 0], Rm[1, 1], t[0], t[1]]
print(f"label pairs {len(pairs)} | used {int(keep.sum())} | residual median {np.median(res[keep]):.2f} m, max {res[keep].max():.2f} m | turn {math.degrees(math.atan2(Rm[1, 0], Rm[0, 0])):.1f} deg")
bad = [(pairs[k][0], round(res[k], 1)) for k in range(len(pairs)) if not keep[k]]
print("label pairs that do not fit (layer label or position differs from the sheet):", bad[:30])
S["geometry"] = [affine_transform(g, Mx) for g in S.geometry]; S = S.set_crs(32644)
S["in_layer"] = S.subdiv_no.isin(set(cur.s))
lumps = cur[cur.kide == SV]; LU = make_valid(unary_union(list(lumps.geometry)))
S["in_lumps_pct"] = [round(100 * g.intersection(LU).area / g.area) for g in S.geometry]
miss = S[(S.subdiv_no != "") & ~S.in_layer]
print(f"register subdivisions on the sheet but not in the layer: {len(miss)} | of them lying >= 80 % inside the unnumbered 124 lumps: {int((miss.in_lumps_pct >= 80).sum())}")
print(miss[["subdiv_no", "in_lumps_pct"]].assign(area=miss.geometry.area.round(1)).sort_values("in_lumps_pct").to_string(index=False))
reg_not_on_sheet = sorted(REG - set(S.subdiv_no) - set(cur.s)); print("register subdivisions neither on the sheet (as read) nor in the layer:", reg_not_on_sheet)
S.to_file("k124_sheet.gpkg", driver="GPKG")
