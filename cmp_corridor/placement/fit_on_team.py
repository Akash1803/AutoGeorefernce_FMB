"""Rigid PDF sheet posed on the TEAM footprint (merged_v1_backup = team CMP layer as delivered) of each survey. -> <vc>_<sv>_team.gpkg"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace"); sys.path.insert(0, ".")
from georef_lib import *
src = open("georef_survey.py", encoding="utf-8").read(); ns = {}; exec(src.split("for sv in sys.argv[2:]:")[0].replace("VC = sys.argv[1]", "VC = None"), ns)
T = gpd.read_file("merged_v1_backup.geojson").to_crs(32644); T["kide"] = T.kide.astype(str); T["sv"] = T.kide.str.split("/").str[0]; T["vc"] = T.village_code.astype(str)
for a in sys.argv[1:]:
    vc, sv = a.split(":")
    c = T[(T.vc == vc) & (T.sv == sv)].copy(); c["s"] = c.subdiv_no.astype(str).str.strip()
    try: p, sk = ns["sheet_plots"](vc, sv)
    except Exception as e: print(vc, sv, "no PDF sheet", e); continue
    if abs(len(p) - len(sk)) > max(2, 0.2 * len(sk)):
        p = sk[["subdiv_no", "geometry"]].copy(); p["subdiv_no"] = p.subdiv_no.astype(str).str.strip().replace({"": "-", "None": "-", "nan": "-"}); frame = "sketch vector"
    else: frame = "PDF"
    Us = make_valid(unary_union(list(p.geometry))); Uc = make_valid(unary_union(list(c.geometry)))
    pairs = [(p[p.subdiv_no == s].geometry.iloc[0].centroid.coords[0], c[c.s == s].geometry.iloc[0].centroid.coords[0]) for s in set(p.subdiv_no) & set(c.s) - {"-", ""} if (p.subdiv_no == s).sum() == 1 and (c.s == s).sum() == 1]
    Mx, rms, rot = ns["fit"](Us, Uc, pairs if len(pairs) >= 3 else None)
    if not (rot <= 20 or rot >= 340):          # sheets here are drawn north-up: keep the turn within +-20 deg
        import math
        from scipy.optimize import minimize
        from scipy.spatial import cKDTree
        A = ns["dens"](Us); cA = np.array(Us.centroid.coords[0]); A0 = A - cA; B = ns["dens"](Uc); tB = cKDTree(B); cB = np.array(Uc.centroid.coords[0])
        def cost(x):
            th, tx, ty = x; ct, st = math.cos(th), math.sin(th); P = A0 @ np.array([[ct, st], [-st, ct]]) + [tx, ty]
            return (np.sqrt((tB.query(P)[0] ** 2).mean()) + np.sqrt((cKDTree(P).query(B)[0] ** 2).mean())) / 2 + 50 * max(0, abs(th) - math.radians(20))
        best = min((minimize(cost, [math.radians(t0), *cB], method="Nelder-Mead", options={"maxiter": 2000}) for t0 in (-15, -5, 0, 5, 15)), key=lambda r: r.fun)
        th, tx, ty = best.x; ct, st = math.cos(th), math.sin(th)
        Mx = [ct, -st, st, ct, tx - (ct * cA[0] - st * cA[1]), ty - (st * cA[0] + ct * cA[1])]; rms = best.fun; rot = math.degrees(th) % 360
        print(f"   {sv}: turn held north-up")
    g = p.copy(); g["geometry"] = [affine_transform(x, Mx) for x in p.geometry]; g = g.set_crs(32644, allow_override=True)
    g["kide"] = [f"{sv}/{s}" if s not in ("-", "") else sv for s in g.subdiv_no]; g["survey_no"] = sv; g["village_code"] = vc; g["poly_id"] = range(1, len(g) + 1)
    g["area_sqm"] = g.geometry.area.round(1); g["area_acre"] = (g.geometry.area / 4046.856).round(4)
    g.to_crs(4326).to_file(f"{vc}_{sv}_team.gpkg", layer=f"{sv}_parcels_auto", driver="GPKG")
    G = make_valid(unary_union(list(g.geometry)))
    print(f"{vc} {sv}: {frame} {len(p)} plots {Us.area:.0f} m2 | team {len(c)} plots {Uc.area:.0f} m2 | rms {rms:.1f} m on {len(pairs)} pairs | turn {rot:.1f} | sheet outside team footprint {G.difference(Uc).area:.0f} m2, team not covered {Uc.difference(G).area:.0f} m2")
