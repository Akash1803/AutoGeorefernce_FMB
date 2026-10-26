"""Gap closer for the re-placed sheets (V2 files): enclosed gaps between placed sheets (30 m morphological closing of the union) are given to a
red sheet, never to the green base or his files. Rules: a gap touching 423 (the tank sheet, drawn 8,100 m2 short of its neighbours) is filled
into 423 as far as 22 m from its outline; any other gap goes to the adjacent red plot with the longest contact. Then T-junctions, despike,
rewrite, and the dual-free topology report. Usage: python gap_close.py"""
import sys, os, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "C:/Users/FAI-AK~1/AppData/Local/Temp/claude/d--Akash-Python/93d4ffec-db1f-4f77-9e3a-dbfb2eece4d9/scratchpad/cmp2")
from georef_lib import *
base = "D:/Projects/CMP_Roads_LPS/FMB_Vector/35_15_028"
BASE = {sv: f"{base}/{sv}_parcels_auto.gpkg" for sv in ("397", "398", "409", "410", "411", "412", "413", "415")}
BASE["414"] = f"{base}/414_parcels_modified.gpkg"
APPROVED = {"416", "418", "419", "420", "423", "142", "141", "140", "424", "425", "451", "188"}  # green since 2026-10-22 evening: never changed again
RED = {f.split("_")[0]: f for f in sorted(os.listdir(".")) if f.endswith("_parcels_auto_V2.gpkg") and f.split("_")[0] not in APPROVED}
ALLV2 = {f.split("_")[0]: f for f in sorted(os.listdir(".")) if f.endswith("_parcels_auto_V2.gpkg")}
G, U = {}, {}
for sv, f in {**BASE, **ALLV2}.items():
    G[sv], U[sv] = load(f)
log = []


def gaps():
    ALL = make_valid(unary_union(list(U.values())))
    closed = ALL.buffer(30, join_style=2).buffer(-30, join_style=2)
    gp = make_valid(closed.difference(ALL))
    out = []
    for g in getattr(gp, "geoms", [gp]):
        if g.area < 20 or g.geom_type != "Polygon":
            continue
        nb = {sv: round(g.buffer(0.1).intersection(U[sv].boundary).length) for sv in U if g.distance(U[sv]) < 0.1}
        nb = {k: v for k, v in nb.items() if v > 2}
        # a real gap is bounded by placed sheets all round; a notch that opens onto land no sheet covers yet (the closing line is its own edge)
        # must stay open - it is another survey's land (e.g. 130 between 139 and 129)
        shared = g.boundary.buffer(0.1).intersection(unary_union([U[sv].boundary for sv in nb])).length if nb else 0.0
        free = 1.0 - min(shared / g.length, 1.0)
        if len(nb) >= 2 and free < 0.10:
            out.append((g, nb))
        elif len(nb) >= 2:
            print(f"notch {g.area:.0f} m2 at ({g.centroid.x:.0f}, {g.centroid.y:.0f}) next to {nb}: {free:.0%} of its edge is open to unplaced land - left open")
    return out


changed = set()
for g, nb in gaps():
    if "423" in nb and "423" in RED:
        fill = make_valid(g.intersection(U["423"].buffer(22, join_style=2)))
        # the pocket at the road junction between 418, 189, 187 and 142 is not the tank: it is 188's land (sheet held), left open
        pocket = unary_union([U[o] for o in ("187", "189") if o in U]).buffer(30, join_style=2)
        fill = make_valid(fill.difference(pocket))
        if fill.is_empty or fill.area < 1:
            log.append(f"gap {g.area:.0f} m2 between {nb}: lies in the 188 pocket, left open")
            continue
        tgt = ("423", G["423"].index[0])
    else:
        best = None
        for sv in nb:
            if sv not in RED:
                continue
            for i, x in zip(G[sv].index, G[sv].geometry):
                L = g.buffer(0.1).intersection(x.boundary).length
                if L > 2 and (best is None or L > best[0]):
                    best = (L, sv, i)
        if best is None:
            log.append(f"gap {g.area:.0f} m2 at ({g.centroid.x:.0f}, {g.centroid.y:.0f}) touches only the fixed base {nb}: left open")
            continue
        fill = g
        tgt = (best[1], best[2])
    sv, i = tgt
    old = G[sv].loc[i, "geometry"]
    new = make_valid(unary_union([old, fill.buffer(0.05, join_style=2)]))
    new = max(getattr(new, "geoms", [new]), key=lambda y: y.area) if new.geom_type != "Polygon" else new
    # never take anything that another sheet, or another plot of the same sheet, already holds
    others = unary_union([U[o] for o in U if o != sv] + [x for j, x in zip(G[sv].index, G[sv].geometry) if j != i])
    new = as_poly(make_valid(new.difference(others)))
    if not new.is_valid:
        new = despike(new)
    G[sv].loc[i, "geometry"] = new
    U[sv] = make_valid(unary_union(list(G[sv].geometry)))
    changed.add(sv)
    if "pdf_check" in G[sv].columns and sv != "423":
        G[sv].loc[i, "pdf_check"] = f"{str(G[sv].loc[i, 'pdf_check'] or '')}; stretched by {fill.area:.0f} m2 onto {', '.join(o for o in nb if o != sv)} to close the gap the sheets leave between them (option A squeeze, the placed side never moves)"
    log.append(f"gap {g.area:.0f} m2 ({2 * g.area / g.length:.1f} m wide) between {nb} -> filled {fill.area:.0f} m2 into {sv}/{G[sv].loc[i, 'subdiv_no']} (plot {old.area:.0f} -> {new.area:.0f} m2)")

# the interior slivers the closing cannot see are the holes of the union: fill those < 60 m2 into the adjacent red plot too
ALL = make_valid(unary_union(list(U.values())))
for pp in getattr(ALL, "geoms", [ALL]):
    for r in pp.interiors:
        h = Polygon(r)
        if h.area < 150:
            best = None
            if "423" in RED and h.distance(U["423"]) < 0.1:
                best = (1e9, "423", G["423"].index[0])
            for sv in (RED if best is None else []):
                for i, x in zip(G[sv].index, G[sv].geometry):
                    L = h.buffer(0.1).intersection(x.boundary).length
                    if L > 0.5 and (best is None or L > best[0]):
                        best = (L, sv, i)
            if best:
                _, sv, i = best
                G[sv].loc[i, "geometry"] = as_poly(make_valid(unary_union([G[sv].loc[i, "geometry"], h])))
                U[sv] = make_valid(unary_union(list(G[sv].geometry)))
                changed.add(sv)
                log.append(f"hole {h.area:.1f} m2 filled into {sv}/{G[sv].loc[i, 'subdiv_no']}")

# T-junctions for the changed sheets against everything near them, then write
for sv in sorted(changed):
    near = [G[o] for o in G if o != sv and U[o].distance(U[sv]) < 3.0]
    G[sv] = tjunctions(G[sv], near)
    dl = []
    G[sv]["geometry"] = [x if x.is_valid else despike(x, dl) for x in G[sv].geometry]
    if dl:
        log.append(f"{sv}: despiked {dl}")
    G[sv]["area_sqm"] = G[sv].geometry.area.round(1)
    G[sv]["area_acre"] = (G[sv].geometry.area / 4046.856).round(4)
    if sv == "423":
        G[sv]["pdf_check"] = "outline from the PDF line work; the drawn outline stops 8,100 m2 short of its imagery-verified neighbours (Puvi's own 423 polygon is 8,755 m2 bigger than the sheet), so the tank polygon was stretched onto them up to 22 m - the neighbours keep their shapes"
    G[sv].to_crs(4326).to_file(RED[sv], layer=f"{sv}_parcels_auto", driver="GPKG")
    U[sv] = make_valid(unary_union(list(G[sv].geometry)))
    log.append(f"{sv}: written ({len(G[sv])} plots, {G[sv].geometry.area.sum() / 4046.856:.3f} ac)")
print("\n".join(log))
seams({sv: G[sv] for sv in G}, "ALL after gap closing")
rem = gaps()
print("remaining enclosed gaps:", [(round(g.area), nb) for g, nb in rem])
