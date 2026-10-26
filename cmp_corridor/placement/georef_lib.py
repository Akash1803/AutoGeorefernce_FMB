"""Shared helpers for the Ponmar batches (PDF line work -> plots, rigid fits, seam warp, checks, attributes, pictures)."""
import math
import numpy as np, geopandas as gpd, pandas as pd, fitz
from shapely.ops import unary_union, polygonize
from shapely.geometry import Point, Polygon, LineString
from shapely.affinity import affine_transform, translate
from shapely import make_valid
from scipy.interpolate import RBFInterpolator
SK = "D:/Projects/CMP_Roads_LPS/FMB_Vector/CMP_FMB_Vector_all_sketches.gpkg"
base = "D:/Projects/CMP_Roads_LPS/FMB_Vector/35_15_028"
PT = 5.669


def kd(g):
    return g.kide.astype(str) if "kide" in g.columns else g.subdiv_no.astype(str)


def rigid(src, dst):
    A = np.array(src)
    B = np.array(dst)
    ca, cb = A.mean(0), B.mean(0)
    H = (A - ca).T @ (B - cb)
    Uu, S, Vt = np.linalg.svd(H)
    R = Vt.T @ Uu.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ Uu.T
    t = cb - R @ ca
    return [R[0, 0], R[0, 1], R[1, 0], R[1, 1], t[0], t[1]], np.hypot(*((R @ A.T).T + t - B).T), math.degrees(math.atan2(R[1, 0], R[0, 0]))


def pt(M, p):
    q = affine_transform(Point(p), M)
    return (q.x, q.y)


def ring(U):
    U = max(getattr(U, "geoms", [U]), key=lambda g: g.area)
    return [tuple(p) for p in np.array(U.exterior.coords)[:-1]]


def load(path):
    g = gpd.read_file(path).to_crs(32644)
    return g, make_valid(unary_union(list(g.geometry)))


def pdf_plots(sv, verbose=True):
    pg = fitz.open(f"D:/Projects/CMP_Roads_LPS/FMB_Sketches/35_15_028/{sv}.pdf")[0]
    W, H = pg.rect.width, pg.rect.height
    lines = []
    for d in pg.get_drawings():
        w = round(d.get("width") or 0, 2)
        col = tuple(d.get("color") or (0, 0, 0))
        if col != (0.0, 0.0, 0.0) or w not in (1.0, 2.0) or str(d.get("dashes")) != "[ 9 0 ] 1":
            continue
        for it in d["items"]:
            if it[0] == "l":
                x1, y1, x2, y2 = it[1].x, it[1].y, it[2].x, it[2].y
                if 0.02 * W < min(x1, x2) and max(x1, x2) < 0.98 * W and 0.03 * H < min(y1, y2) and max(y1, y2) < 0.97 * H:
                    lines.append(LineString([(x1 / PT, -y1 / PT), (x2 / PT, -y2 / PT)]))
    polys = [p for p in polygonize(unary_union(lines)) if p.area > 1.0]
    sk = gpd.read_file(SK, layer=f"35_15_028_{sv}")
    sk["subdiv_no"] = sk.subdiv_no.astype(str).str.strip()
    big, refp = max(polys, key=lambda p: p.area), max(sk.geometry, key=lambda g: g.area)
    polys = [translate(p, refp.centroid.x - big.centroid.x, refp.centroid.y - big.centroid.y) for p in polys]
    rows = []
    for p in polys:
        hit = [r.subdiv_no for _, r in sk.iterrows() if r.geometry.intersection(p).area > 0.6 * p.area and r.subdiv_no]
        rows.append({"subdiv_no": hit[0] if hit else "-", "geometry": p})
    g = gpd.GeoDataFrame(rows, geometry="geometry")
    if verbose:
        print(f"{sv}: {len(g)} pdf plots (vector {len(sk)}), unnamed {int((g.subdiv_no == '-').sum())}, area {unary_union(list(g.geometry)).area:.0f} vs vector {unary_union(list(sk.geometry)).area:.0f} m2")
    return g, sk


def edge_vertex(r, length, which="start", tol=1.0, near=None):
    """vertex of the ring edge with the given length (start or end, in ring order); `near` picks among several candidates"""
    n = len(r)
    ks = [k for k in range(n) if abs(math.dist(r[k], r[(k + 1) % n]) - length) < tol]
    assert ks, f"no edge of {length} m"
    if near is not None and len(ks) > 1:
        ks.sort(key=lambda k: math.dist(r[k], near))
    k = ks[0]
    return r[k] if which == "start" else r[(k + 1) % n]


def path_along(r, a, b, direction=1):
    """ring path from a to b stepping `direction`"""
    i, j = r.index(a), r.index(b)
    path = [r[i]]
    k = i
    while k != j:
        k = (k + direction) % len(r)
        path.append(r[k])
    return path


def shorter_path(r, a, b):
    p1, p2 = path_along(r, a, b, 1), path_along(r, a, b, -1)
    return min((p1, p2), key=lambda p: LineString(p).length)


def point_along(path, dist):
    L = LineString(path)
    q = L.interpolate(min(dist, L.length))
    return (q.x, q.y)


def snap_to(gdf, paths, fixed_pts, tol=0.5):
    Ls = [LineString(p) for p in paths]
    moved = []
    Uo = make_valid(unary_union([make_valid(x) for x in gdf.geometry]))
    outline = unary_union([pp.exterior for pp in getattr(Uo, "geoms", [Uo]) if pp.geom_type == "Polygon"])

    def fx(poly):
        poly = as_poly(poly)
        c = list(poly.exterior.coords)[:-1]
        out = []
        for (x, y) in c:
            p = Point(x, y)
            best = min(Ls, key=lambda L: L.distance(p))
            d = best.distance(p)
            if 0 < d < tol and outline.distance(p) < 0.01:
                q = best.interpolate(best.project(p))
                out.append((q.x, q.y))
                moved.append(round(d, 3))
            else:
                out.append((x, y))
        res = []
        for i in range(len(out)):
            a, b = out[i], out[(i + 1) % len(out)]
            res.append(a)
            seg = LineString([a, b])
            inner = []
            for L, path in zip(Ls, paths):
                if L.distance(Point(a)) < 0.002 and L.distance(Point(b)) < 0.002:
                    ta, tb = L.project(Point(a)), L.project(Point(b))
                    inner += [v for v in path if min(ta, tb) + 0.01 < L.project(Point(v)) < max(ta, tb) - 0.01]
            inner += [v for v in fixed_pts if seg.distance(Point(v)) < 0.003 and math.dist(v, a) > 0.005 and math.dist(v, b) > 0.005 and v not in inner]
            inner.sort(key=lambda v: math.dist(v, a))
            res.extend(inner)
        return Polygon(res + [res[0]], [list(r.coords) for r in poly.interiors])
    g = gdf.copy()
    g["geometry"] = [fx(x) for x in gdf.geometry]
    return g, moved


def warp_snap(gdf, paths, tol):
    """outline vertices within tol of a path move onto it (corner to corner when a path corner is within reach); every other vertex follows smoothly; then exact snap + T-junctions"""
    Ls = [LineString(p) for p in paths]
    Uo = make_valid(unary_union([make_valid(x) for x in gdf.geometry]))
    outline_pts = []
    for pp in getattr(Uo, "geoms", [Uo]):
        outline_pts += [tuple(v) for v in pp.exterior.coords[:-1]]
    path_vertices = [tuple(v) for p in paths for v in p]
    src, dst, moved = [], [], []
    for v in outline_pts:
        p = Point(v)
        vn = min(path_vertices, key=lambda w: math.dist(w, v))
        best = min(Ls, key=lambda L: L.distance(p))
        d = best.distance(p)
        if 0 < math.dist(vn, v) < tol:
            src.append(v)
            dst.append(vn)
            moved.append(round(math.dist(vn, v), 3))
        elif 0 < d < tol:
            q = best.interpolate(best.project(p))
            src.append(v)
            dst.append((q.x, q.y))
            moved.append(round(d, 3))
        else:
            src.append(v)
            dst.append(v)
    tps = RBFInterpolator(np.array(src), np.array(dst), kernel="thin_plate_spline", smoothing=0.0)

    def wp(poly):
        poly = as_poly(poly)
        ext = tps(np.array(poly.exterior.coords)[:-1])
        ints = [tps(np.array(r.coords)[:-1]) for r in poly.interiors]
        return Polygon([tuple(x) for x in ext] + [tuple(ext[0])], [[tuple(x) for x in r] + [tuple(r[0])] for r in ints])
    g = gdf.copy()
    g["geometry"] = [wp(x) for x in gdf.geometry]
    g["geometry"] = [x if x.is_valid else despike(x) for x in g.geometry]
    g, _ = snap_to(g, paths, [], tol=0.02)
    g["geometry"] = [x if x.is_valid else despike(x) for x in g.geometry]
    return g, moved


def despike(poly, log=None):
    """remove vertices where the ring doubles back on itself (turn within 2 deg of a reversal) until the polygon is valid; last resort make_valid (largest part)"""
    if poly.is_valid:
        return poly
    c = list(poly.exterior.coords)[:-1]
    changed = True
    while changed and not Polygon(c + [c[0]], [list(r.coords) for r in poly.interiors]).is_valid and len(c) > 3:
        changed = False
        for i in range(len(c)):
            a, v, b = c[i - 1], c[i], c[(i + 1) % len(c)]
            d1 = math.atan2(a[1] - v[1], a[0] - v[0])
            d2 = math.atan2(b[1] - v[1], b[0] - v[0])
            turn = abs((math.degrees(d1 - d2) + 180) % 360 - 180)
            if turn < 2.0 or math.dist(a, v) < 0.002 or math.dist(v, b) < 0.002:
                c.pop(i)
                changed = True
                if log is not None:
                    log.append(("despike", tuple(np.round(v, 2))))
                break
    out = Polygon(c + [c[0]], [list(r.coords) for r in poly.interiors])
    if not out.is_valid:
        fixed = make_valid(out)
        out = max(getattr(fixed, "geoms", [fixed]), key=lambda g: g.area if g.geom_type == "Polygon" else 0)
        if log is not None:
            log.append(("make_valid", round(out.area, 1)))
    return out


def as_poly(x):
    """largest polygon of a multi-part result (overlay leftovers)"""
    if x.geom_type == "Polygon":
        return x
    parts = [g for g in getattr(x, "geoms", [x]) if g.geom_type == "Polygon"]
    return max(parts, key=lambda g: g.area) if parts else x


def tjunctions(gdf, others):
    pts = [p for g in others for geom in g.geometry for p in as_poly(geom).exterior.coords[:-1]]

    def fx(poly):
        poly = as_poly(poly)
        c = list(poly.exterior.coords)[:-1]
        res = []
        for i in range(len(c)):
            a, b = c[i], c[(i + 1) % len(c)]
            res.append(a)
            seg = LineString([a, b])
            inner = [v for v in pts if seg.distance(Point(v)) < 0.003 and math.dist(v, a) > 0.005 and math.dist(v, b) > 0.005]
            inner.sort(key=lambda v: math.dist(v, a))
            res.extend(inner)
        out = Polygon(res + [res[0]], [list(r.coords) for r in poly.interiors])
        return out if out.is_valid else despike(out)
    g = gdf.copy()
    g["geometry"] = [fx(x) for x in gdf.geometry]
    return g


def seams(gdfs, label):
    A = pd.concat([g.assign(lab=n + ":" + kd(g))[["lab", "geometry"]] for n, g in gdfs.items()], ignore_index=True)
    bad = A.lab[~A.is_valid].tolist()
    A["geometry"] = [make_valid(x) for x in A.geometry]
    sidx = A.sindex
    ov = [(A.lab[i], A.lab[j], round(g.intersection(A.geometry[j]).area, 3)) for i, g in enumerate(A.geometry) for j in sidx.query(g, predicate="intersects") if j > i and g.intersection(A.geometry[j]).area > 0.01]
    U = make_valid(unary_union(list(A.geometry)))
    holes = [(round(Polygon(r).area, 3), tuple(np.round(Polygon(r).centroid.coords[0], 1)), [A.lab[i] for i, g in enumerate(A.geometry) if g.distance(Polygon(r)) < 0.001]) for p in getattr(U, "geoms", [U]) for r in p.interiors if Polygon(r).area > 0.01]
    G = U.buffer(3, join_style=2).buffer(-3, join_style=2).difference(U)
    sl = sorted([(round(x.area, 2), tuple(np.round(x.centroid.coords[0], 1))) for x in getattr(G, "geoms", [G]) if x.area > 0.3], key=lambda t: -t[0])
    print(f"{label}: invalid {bad} | overlaps {sorted(ov, key=lambda t: -t[2])[:8]} | holes {holes[:6]} | slivers {sl[:6]}")
    return ov, holes, sl


def area_report(g, sk, sv):
    sa = {r.subdiv_no: r.geometry.area for _, r in sk.iterrows() if r.subdiv_no}
    dev = sorted([(r.subdiv_no, round((r.geometry.area / sa[r.subdiv_no] - 1) * 100, 1)) for _, r in g.iterrows() if r.subdiv_no in sa], key=lambda t: abs(t[1]))
    print(f"{sv}: plot area vs sheet, largest changes {dev[-4:]}")


def finish(g, sv, placed_by, pdf_check="number matches the PDF; shape from the PDF line work"):
    tmpl = gpd.read_file(f"{base}/418_parcels_auto.gpkg").iloc[0]
    g = g.copy()
    g["subdiv_no"] = [s if s else "-" for s in g.subdiv_no]
    g["survey_no"] = sv
    g["kide"] = [f"{sv}/" + s if s != "-" else sv for s in g.subdiv_no]
    g["area_sqm"] = g.geometry.area.round(1)
    g["area_acre"] = (g.geometry.area / 4046.856).round(4)
    for c in ("village", "village_code", "district", "taluk", "sheet_scale"):
        g[c] = tmpl[c]
    g["giscode"] = str(tmpl["giscode"])[:-3] + sv if isinstance(tmpl["giscode"], str) else tmpl["giscode"]
    g["source_pdf"] = f"FMB_Sketches\\35_15_028\\{sv}.pdf"
    g["poly_id"] = range(1, len(g) + 1)
    g["unique_id"] = [f"35_15_028_{sv}/" + s for s in g.subdiv_no]
    g["pdf_check"] = pdf_check
    g["placed_by"] = placed_by
    g["status"] = "red - for review"
    front = ["kide", "survey_no", "subdiv_no", "area_acre", "area_sqm", "village", "village_code", "district", "taluk", "giscode", "sheet_scale", "source_pdf", "poly_id", "unique_id", "pdf_check", "placed_by", "status"]
    return g[front + ["geometry"]].set_crs(32644, allow_override=True)


def picture(name, title, layers, sat_file="sat_south2.png", X0=409900.0, Y1=1421750.0, Rr=0.25, pad=30):
    from PIL import Image, ImageDraw, ImageFont
    Image.MAX_IMAGE_PIXELS = None
    sat = Image.open(sat_file).convert("RGB")
    U = unary_union([make_valid(unary_union(list(g.geometry))) for g, col in layers if col != "ctx"])
    ext = U.bounds
    cx, cy, half = (ext[0] + ext[2]) / 2, (ext[1] + ext[3]) / 2, max(ext[2] - ext[0], ext[3] - ext[1]) / 2 + pad
    img = sat.crop((int((cx - half - X0) / Rr), int((Y1 - (cy + half)) / Rr), int((cx + half - X0) / Rr), int((Y1 - (cy - half)) / Rr)))
    dr = ImageDraw.Draw(img)
    font = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 14)
    small = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 10)
    ox, oy = (cx - half - X0) / Rr, (Y1 - (cy + half)) / Rr

    def q(x, y, *_):
        return ((x - X0) / Rr - ox), ((Y1 - y) / Rr - oy)
    for g, col in layers:
        c_ = (255, 120, 0) if col == "ctx" else col
        for k, geom in zip(kd(g), g.geometry):
            for poly in getattr(geom, "geoms", [geom]):
                dr.line([q(*p) for p in poly.exterior.coords], fill=c_, width=2)
            if geom.area > 150:
                c = geom.representative_point()
                dr.text(q(c.x, c.y), str(k).split("/")[-1], fill=c_, font=small)
    dr.text((6, 6), title, fill=(255, 255, 255), font=font)
    img.save(name)
    print(name, img.size, "extent", [round(v) for v in ext])
