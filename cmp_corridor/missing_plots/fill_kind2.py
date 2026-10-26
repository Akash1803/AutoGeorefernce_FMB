"""Add PDF plots that the vector conversion left out (kind 2) into CMP_two_Parts: the warped sheet plot (k2_warped.gpkg, label-anchored
pose) is clipped to the ONE host parcel holding >= 85 % of it; the host loses exactly that area; nothing else changes.
Approved list only (Akash 2026-10-26). -> fix_two_fill.json (host change + new plot) + report"""
import sys, json, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd
from shapely import make_valid
from shapely.geometry import Polygon
norm = lambda s: re.sub(r"[^0-9A-Z]", "", str(s).upper())
APPROVED = [("35_11_047", "344", "2"), ("35_11_047", "368", "2A6"), ("35_11_047", "371", "3F"), ("35_11_047", "767", "4"), ("35_11_047", "898", "1B1"),
            ("35_11_047", "905", "1A"), ("35_11_047", "943", "4B1A4"), ("35_11_101", "255", "10A"), ("35_11_117", "36", "3B2"), ("35_11_117", "38", "3A1A1A1")]
W = gpd.read_file("k2_warped.gpkg")
A = gpd.read_file(r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"); A["vc"] = A.village_code.astype(str); A["kide"] = A.kide.astype(str)
G = {i: make_valid(g) for i, g in zip(A.index, A.geometry)}
out = {"change": [], "add": [], "delete": []}; rep = []
for vc, sv, num in APPROVED:
    w = W[(W.village == vc) & (W.survey == sv)]; hit = w[[norm(num) in {norm(t.rstrip("*")) for t in n.split("|") if t} for n in w.numbers]]
    if not len(hit): rep.append((sv, num, "sheet plot not found")); continue
    g = make_valid(hit.geometry.iloc[0])
    hosts = sorted([(G[i].intersection(g).area / g.area, i) for i in A.index[A.intersects(g)]], reverse=True)
    if not hosts or hosts[0][0] < 0.85: rep.append((sv, num, f"not inside one parcel ({hosts[0][0]:.0%})")); continue
    hi = hosts[0][1]; host = G[hi]
    new = make_valid(g.intersection(host)); new = max([z for z in getattr(new, "geoms", [new]) if z.geom_type == "Polygon"], key=lambda z: z.area)
    rest = make_valid(host.difference(new.buffer(0.001, join_style=2)))
    parts = sorted([z for z in getattr(rest, "geoms", [rest]) if z.geom_type == "Polygon"], key=lambda z: -z.area)
    if len(parts) != 1 and sum(z.area for z in parts[1:]) > 1.0: rep.append((sv, num, f"would split {A.kide[hi]} into {len(parts)} parts")); continue
    rest = Polygon(parts[0].exterior, [r for r in parts[0].interiors if Polygon(r).area > 0.5])
    G[hi] = rest
    out["change"].append({"vc": vc, "old_kide": A.kide[hi], "old_pt": list(host.representative_point().coords[0]), "wkt": rest.wkt,
                          "attrs": {"area_sqm": round(rest.area, 1), "area_acre": round(rest.area / 4046.856, 4), "merge_note": (str(A.merge_note[hi] or "") + f"; 2026-10-26: plot {sv}/{num} (missing in the vector conversion) cut out of it from the PDF").strip("; ")}})
    at = {k: (None if str(v) in ("nan", "None", "NaT") else v) for k, v in A.loc[hi].drop(labels=["geometry", "vc"]).to_dict().items()}
    for k in list(at):
        if k.startswith("tamilnilam_") or k.startswith("reginet_"): at[k] = None
    at.update({"kide": f"{sv}/{num}", "survey_no": sv, "subdiv_no": num, "unique_id": f"{vc}_{sv}/{num}", "area_sqm": round(new.area, 1), "area_acre": round(new.area / 4046.856, 4),
               "source": f"georeferenced from the {sv} PDF (plot missing in the vector conversion)", "status": "red - for review", "review": "rigid placement - for review",
               "pdf_check": f"plot {num} printed on the {sv} PDF; sheet posed on the numbered plots around it; cut from {A.kide[hi]}", "merge_note": "", "placed_by": "auto 2026-10-26"})
    out["add"].append({"wkt": new.wkt, "attrs": at})
    rep.append((sv, num, f"added {new.area:.0f} m2 (PDF {g.area:.0f}), cut from {A.kide[hi]} {host.area:.0f} -> {rest.area:.0f}"))
for r in rep: print(r)
json.dump(out, open("fix_two_fill.json", "w"), default=str); print("changes", len(out["change"]), "adds", len(out["add"]))
