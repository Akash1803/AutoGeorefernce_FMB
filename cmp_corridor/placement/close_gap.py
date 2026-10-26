"""Close the gaps around one georeferenced sheet by RIGID placement (no new shapes, no reshaping):
1. the full PDF sheet (all plots, true scale) is posed onto the sheet's present rows in CMP_my_Part (rms printed);
2. icp_lib.close_onto slides/turns it (turn <= CLAMP deg, move <= MAXMOVE m) against the outlines of the parcels around it;
3. reports gap and overlap around the sheet before/after, writes <vc>_<sv>_closed.gpkg (for apply_replacement.py) and
   closed_<sv>.png (Google Satellite: orange = neighbours, red = sheet now, cyan = closed position).
Usage: python close_gap.py <vc> <sv> [maxmove=15] [clamp=3] [ignore_sv,...]"""
import sys, math, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "."); sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
from georef_lib import *
from icp_lib import close_onto
from shapely.geometry import box
from PIL import Image, ImageDraw, ImageFont
from autogeoref import satcheck
vc, sv = sys.argv[1], sys.argv[2]
MAXMOVE = float(sys.argv[3]) if len(sys.argv) > 3 else 15.0
CLAMP = float(sys.argv[4]) if len(sys.argv) > 4 else 3.0
IGN = set(sys.argv[5].split(",")) if len(sys.argv) > 5 and sys.argv[5] else set()
src = open("georef_survey.py", encoding="utf-8").read()
ns = {}; exec(src.split("for sv in sys.argv[2:]:")[0].replace("VC = sys.argv[1]", "VC = None"), ns)
M = ns["M"]; M["vc"] = M.village_code; M["sv"] = M.survey
cur = M[(M.vc == vc) & (M.sv == sv)]
assert len(cur) and cur.source.astype(str).str.startswith("georef").all(), f"{sv}: not a georeferenced sheet in CMP_my_Part"
Uc = make_valid(unary_union(list(cur.geometry)))
# 1. full sheet posed onto the present rows
try:
    p, sk = ns["sheet_plots"](vc, sv)
    if abs(len(p) - len(cur)) > max(2, 0.2 * len(cur)): raise ValueError(f"PDF gives {len(p)} plots, layer {len(cur)}")
    Us = make_valid(unary_union(list(p.geometry)))
    c2 = cur.copy(); c2["s"] = c2.subdiv_no.astype(str).str.strip()
    pairs = [(p[p.subdiv_no == s].geometry.iloc[0].centroid.coords[0], c2[c2.s == s].geometry.iloc[0].centroid.coords[0])
             for s in set(p.subdiv_no) & set(c2.s) - {"-", ""} if (p.subdiv_no == s).sum() == 1 and (c2.s == s).sum() == 1]
    Mx, rms, rot = ns["fit"](Us, Uc, pairs if len(pairs) >= 3 else None)
    g = p.copy(); g["geometry"] = [affine_transform(x, Mx) for x in p.geometry]; g = g.set_crs(32644, allow_override=True)
    frame = f"PDF sheet posed on the present rows (rms {rms:.1f} m)"
    if rms > 3.0: raise ValueError(f"pose rms {rms:.1f} m")
except Exception as e:
    g = cur[["subdiv_no", "geometry"]].copy().reset_index(drop=True); frame = f"present rows as they are ({e})"
g["subdiv_no"] = g.subdiv_no.astype(str).str.strip().replace({"": "-", "None": "-", "nan": "-"})
G0 = make_valid(unary_union(list(g.geometry)))
# 2. close onto the neighbours
nbh = M[M.intersects(G0.buffer(40)) & ~((M.vc == vc) & (M.sv == sv)) & ~((M.vc == vc) & M.sv.isin(IGN)) & (M.geometry.area < 300000)]
nb = [make_valid(x) for x in nbh.geometry]
g1, info = close_onto(g, nb, clamp=CLAMP, maxmove=MAXMOVE)
G1 = make_valid(unary_union(list(g1.geometry)))
geo0 = make_valid(unary_union([x for x, s in zip(nbh.geometry, nbh.source.astype(str)) if s.startswith("georef")])) if len(nbh) else Polygon()
cmp0 = make_valid(unary_union([x for x, s in zip(nbh.geometry, nbh.source.astype(str)) if not s.startswith("georef")])) if len(nbh) else Polygon()
NBU = make_valid(unary_union(nb)) if nb else Polygon()
def score(U):   # open ground left next to the sheet + overlap (placed sheets count 3x: the new sheet would be cut there)
    A = NBU.union(U); cl = A.buffer(8, join_style=2).buffer(-8, join_style=2); d = make_valid(cl.difference(A)).intersection(U.buffer(12))
    return d.area + 3 * U.intersection(geo0).area + U.intersection(cmp0).area
import itertools
from shapely.affinity import rotate as _rot, translate as _tr
if "move_m" in info or info.get("note") == "no neighbours":
    base = G1; best = (score(base), 0, 0, 0)
    for step, span, astep in ((1.0, 4, 0.5), (0.25, 1, 0.25)):
        b0 = best
        for dx, dy, da in itertools.product(np.arange(-span, span + 1e-9, step), np.arange(-span, span + 1e-9, step), np.arange(-2 * astep, 2 * astep + 1e-9, astep)):
            X, Y, A_ = b0[1] + dx, b0[2] + dy, b0[3] + da
            sc = score(_tr(_rot(base, A_, origin=base.centroid), X, Y))
            if sc < best[0] - 0.5: best = (sc, X, Y, A_)
    _, X, Y, A_ = best; c0 = base.centroid
    g1 = g1.copy(); g1["geometry"] = [_tr(_rot(x, A_, origin=c0), X, Y) for x in g1.geometry]; G1 = make_valid(unary_union(list(g1.geometry)))
    info["refine"] = f"dx {X:+.2f} dy {Y:+.2f} turn {A_:+.2f}"
    tot = G1.centroid; info["total_move_m"] = round(math.hypot(tot.x - G0.centroid.x, tot.y - G0.centroid.y), 1)
    if info["total_move_m"] > MAXMOVE: print("REJECT: total move", info["total_move_m"]); g1 = g; G1 = G0
geo = make_valid(unary_union([x for x, s in zip(nbh.geometry, nbh.source.astype(str)) if s.startswith("georef")])) if len(nbh) else Polygon()
cmp_ = make_valid(unary_union([x for x, s in zip(nbh.geometry, nbh.source.astype(str)) if not s.startswith("georef")])) if len(nbh) else Polygon()
def gaps(U):
    A = make_valid(unary_union(nb + [U])); cl = A.buffer(8, join_style=2).buffer(-8, join_style=2)
    d = make_valid(cl.difference(A)); win = U.buffer(12)
    return sum(q.area for q in getattr(d, "geoms", [d]) if q.geom_type == "Polygon" and q.area > 5 and q.intersects(win))
for nm, U in (("now", Uc), ("full sheet at present pose", G0), ("closed", G1)):
    print(f"  {nm:28s}: gap around {gaps(U):7.0f} m2 | overlap georef sheets {U.intersection(geo).area:6.0f} m2 | overlap CMP parcels {U.intersection(cmp_).area:6.0f} m2 | area {U.area:.0f}")
print(f"{vc} {sv}: {frame}; close: {info}")
g1["kide"] = [f"{sv}/{s}" if s not in ("-", "") else sv for s in g1.subdiv_no]; g1["survey_no"] = sv; g1["village_code"] = vc
g1["poly_id"] = range(1, len(g1) + 1); g1["area_sqm"] = g1.geometry.area.round(1); g1["area_acre"] = (g1.geometry.area / 4046.856).round(4)
g1.to_crs(4326).to_file(f"{vc}_{sv}_closed.gpkg", layer=f"{sv}_parcels_auto", driver="GPKG")
# 3. picture
x0, y0, x1, y1 = G1.union(Uc).buffer(35).bounds
jpg = os.path.abspath(f"satc_{vc}_{sv}.jpg")
if not os.path.exists(jpg): satcheck.render_from_qgis((x0, y0, x1, y1), jpg)
im = Image.open(jpg).convert("RGB"); W, H = im.size; jx0, jy1 = x0, y1; rx = (x1 - x0) / W
k = max(1, int(1600 / max(W, H))) if max(W, H) < 1600 else 1
im = im.resize((W * k, H * k)); dr = ImageDraw.Draw(im); q = lambda x, y, *_: ((x - jx0) / rx * k, (jy1 - y) / rx * k)
for x in nbh.geometry:
    for pp in getattr(x, "geoms", [x]):
        if pp.geom_type == "Polygon": dr.line([q(*c) for c in pp.exterior.coords], fill=(255, 150, 0), width=1)
for x in cur.geometry:
    for pp in getattr(x, "geoms", [x]): dr.line([q(*c) for c in pp.exterior.coords], fill=(255, 0, 0), width=2)
def gap_polys(U):
    A = make_valid(unary_union(nb + [U])); cl = A.buffer(8, join_style=2).buffer(-8, join_style=2); d = make_valid(cl.difference(A))
    return [q_ for q_ in getattr(d, "geoms", [d]) if q_.geom_type == "Polygon" and q_.area > 5]
for gp in gap_polys(Uc): dr.polygon([q(*c) for c in gp.exterior.coords], outline=(255, 0, 255))
for gp in gap_polys(G1): dr.line([q(*c) for c in gp.exterior.coords], fill=(255, 255, 255), width=1)
f = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 13)
for x, l in zip(g1.geometry, g1.subdiv_no):
    dr.line([q(*c) for c in x.exterior.coords], fill=(0, 255, 255), width=2)
    if x.area > 150: cc = x.representative_point(); dr.text(q(cc.x, cc.y), l, fill=(255, 255, 0), font=f)
dr.text((6, 6), f"{sv}: orange = neighbours, red = now (magenta = its gaps), cyan = closed (white = gaps left) ({info.get('move_m', '-')} m, {info.get('turn_deg', '-')} deg)", fill=(255, 255, 255), font=ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 16))
im.save(f"closed_{sv}.png"); print(f"  -> {vc}_{sv}_closed.gpkg, closed_{sv}.png")
