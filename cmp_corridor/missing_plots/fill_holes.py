"""Close the remaining holes of CMP_two_Parts that are missing PDF plots/faces (vector-conversion losses): each hole is cut along the
warped PDF plot boundaries (label-anchored warp, chk/k2_warped.gpkg) and every piece becomes a parcel of that survey - numbered only when
the PDF number is certain, in the Tamilnilam register and not yet in the layer; otherwise an unnumbered piece with the read in pdf_check.
Nothing existing changes. Usage: fill_holes.py <gaps.csv> <row:survey,row:survey,...>  -> fix_two_holes.json"""
import sys, json, re, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import geopandas as gpd, pandas as pd, numpy as np
from shapely import make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union, polygonize
F = r"D:\Projects\CMP_Roads_LPS\FMB_Georeferenced\CMP_two_Parts.geojson"
norm = lambda s: re.sub(r"[^0-9A-Z]", "", str(s).upper())
P = gpd.read_file(F); P = P[P.geometry.notna()].reset_index(drop=True); P["kide"] = P.kide.astype(str); P["vc"] = P.village_code.astype(str)
G = [make_valid(g) for g in P.geometry]; U = unary_union(G)
holes = [Polygon(r) for pp in getattr(U, "geoms", [U]) if pp.geom_type == "Polygon" for r in pp.interiors]
D = pd.read_csv(sys.argv[1])
W = pd.concat([gpd.read_file("chk_warped.gpkg"), gpd.read_file("k2_warped.gpkg")], ignore_index=True).drop_duplicates(subset=["village", "survey", "pdf_plot"])
T = pd.read_excel("D:/Projects/CMP_Roads_LPS/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
REG = T.groupby([T.vc, T.fmb_survey_no.map(norm)]).sub_division.apply(lambda x: {norm(v) for v in x} - {""}).to_dict()
out = {"change": [], "add": [], "delete": []}; rep = []
for item in sys.argv[2].split(","):
    r, sv = item.split(":"); row = D.iloc[int(r)]
    h = min(holes, key=lambda z: abs(z.area - row.area_m2) + z.distance(Polygon([(row.x, row.y)] * 4).buffer(0.1)))
    if os.environ.get("FH_CLIP"):      # only the part of the hole inside the corridor buffer (+5 m) is filled
        BUF = unary_union(gpd.read_file(r"D:\Projects\CMP_Roads_LPS\CMP_Roads_Buffer.gpkg", layer="road_buffer_30m").to_crs(32644).geometry.values).buffer(float(os.environ["FH_CLIP"]), join_style=2)
        hc = make_valid(h.intersection(BUF)); h = max([z for z in getattr(hc, "geoms", [hc]) if z.geom_type == "Polygon"], key=lambda z: z.area)
        print(f"hole {r}: clipped to the corridor -> {h.area:.0f} m2")
    nb = P[P.intersects(h.buffer(0.02)) & (P.kide.str.split("/").str[0] == sv)]
    assert len(nb), f"hole {r}: no {sv} parcel beside it"
    vc = nb.vc.iloc[0]; have = set(P.kide[P.vc == vc].map(lambda k: norm(k.split("/")[1]) if "/" in k else ""))
    w = W[(W.village == vc) & (W.survey == sv)]; w = w[[g.intersection(h).area > 0.05 * h.area for g in w.geometry]]
    lines = [h.exterior] + [g.boundary for g in w.geometry]
    faces = [f for f in polygonize(unary_union(lines)) if f.representative_point().within(h) and f.area > 0.01]
    # face -> owning PDF plot (biggest share); faces without a plot join the neighbouring face with the longest shared edge
    own = {}
    for k, f in enumerate(faces):
        best = max(((g.intersection(f).area / f.area, p) for p, g in zip(w.pdf_plot, w.geometry)), default=(0, None))
        own[k] = best[1] if best[0] > 0.5 else None
    groups = {}
    for k, f in enumerate(faces):
        if own[k] is None:
            cands = [(f.buffer(0.01).intersection(faces[j]).area, j) for j in range(len(faces)) if j != k and own[j] is not None]
            own[k] = own[max(cands)[1]] if cands else "-"
        groups.setdefault(own[k], []).append(f)
    # unnumbered pieces under 80 m2 join the piece of the same hole they share the longest edge with (no slivers)
    def numbered(p):
        if p == "-": return False
        nums = [t for t in str(w.numbers[w.pdf_plot == p].iloc[0]).split("|") if t]; cert = [norm(t) for t in nums if not t.endswith("*")]
        return len(nums) == 1 and len(cert) == 1 and cert[0] in REG.get((vc, norm(sv)), set()) and cert[0] not in have
    reads = {}
    changed = True
    while changed and len(groups) > 1:
        changed = False
        for p in list(groups):
            gp = unary_union(groups[p])
            if gp.area < 80 and not numbered(p):
                others = [(gp.buffer(0.01).intersection(unary_union(groups[q])).area, q) for q in groups if q != p]
                q = max(others)[1]; groups[q] += groups.pop(p); reads.setdefault(q, []).append(p); changed = True; break
    tmpl = nb.iloc[0]
    for p, fs in groups.items():
        g = make_valid(unary_union(fs)).buffer(0)
        g = max([z for z in getattr(g, "geoms", [g]) if z.geom_type == "Polygon"], key=lambda z: z.area) if g.geom_type != "Polygon" else g
        nums = [] if p == "-" else [t for t in str(w.numbers[w.pdf_plot == p].iloc[0]).split("|") if t]
        certain = [norm(t) for t in nums if not t.endswith("*")]
        reg = REG.get((vc, norm(sv)), set())
        num = certain[0] if len(nums) == 1 and len(certain) == 1 and certain[0] in reg and certain[0] not in have else ""
        if num: have.add(num)
        kide = f"{sv}/{num}" if num else sv
        read = ("no number printed (layout face)" if not nums else "printed number read as " + ", ".join(nums).replace("*", " (uncertain)"))
        extra = "".join(f"; also PDF plot {q} ({str(w.numbers[w.pdf_plot == q].iloc[0]) or 'no number'})" for q in reads.get(p, []))
        note = f"PDF plot {p} of sheet {sv} missing in the vector conversion; {read}{extra}; " + (f"numbered {num} - in the register, not yet in the layer" if num else "left unnumbered (number uncertain, already used, or not printed)") + " - added 2026-10-26, cut to the gap"
        at = {c: (None if str(tmpl[c]) in ("nan", "None", "NaT") else tmpl[c]) for c in P.columns if c not in ("geometry", "vc")}
        for c in list(at):
            if c.startswith("tamilnilam_") or c.startswith("reginet_"): at[c] = None
        at.update({"kide": kide, "survey_no": sv, "subdiv_no": num, "unique_id": f"{vc}_{kide}" + ("" if num else "/-"), "area_sqm": round(g.area, 1), "area_acre": round(g.area / 4046.856, 4),
                   "source": f"georeferenced from the {sv} PDF (plot missing in the vector conversion)", "status": "red - for review", "review": "PDF face placed by the label-anchored warp, cut to the gap - for review",
                   "pdf_check": note, "merge_note": "", "placed_by": "auto 2026-10-26 (PDF face clipped to the gap)", "poly_id": None})
        out["add"].append({"wkt": g.wkt, "attrs": at}); rep.append((r, sv, p, kide, round(g.area, 1), read))
    covered = unary_union([f for fs in groups.values() for f in fs]).area
    rep.append((r, sv, "hole %.1f m2 -> covered %.1f m2" % (h.area, covered), "", "", ""))
for x in rep: print(x)
json.dump(out, open("fix_two_holes.json", "w"), default=str); print("adds", len(out["add"]))
