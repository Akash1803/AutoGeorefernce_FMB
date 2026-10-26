"""Rigid closing of one sheet onto already placed neighbours (the chain_v2 ICP, reusable): facing-edge matches only, rotation clamped to
+-CLAMP deg of the start pose, move capped. close_onto(gdf, neighbour_polys, clamp=3.0, maxmove=25.0) -> (gdf_moved, info)"""
import math
import numpy as np
from shapely.affinity import affine_transform
from shapely.geometry.polygon import orient
from shapely.ops import unary_union
from shapely import make_valid
from scipy.spatial import cKDTree


def dens_pts(geom, step=2.0):
    pts, nrm = [], []
    for pp in getattr(geom, "geoms", [geom]):
        if pp.geom_type != "Polygon":
            continue
        L = orient(pp, 1.0).exterior
        n = max(int(L.length / step), 8)
        for i in range(n):
            a = L.interpolate(i / n, normalized=True).coords[0]
            b = L.interpolate((i + 0.5) / n, normalized=True).coords[0]
            c = L.interpolate(((i - 0.5) % n) / n, normalized=True).coords[0]
            dx, dy = b[0] - c[0], b[1] - c[1]
            h = math.hypot(dx, dy) or 1.0
            pts.append(a)
            nrm.append((dy / h, -dx / h))
    return np.array(pts), np.array(nrm)


def kabsch(src, dst):
    A = np.array(src); B = np.array(dst)
    ca, cb = A.mean(0), B.mean(0)
    H = (A - ca).T @ (B - cb)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T
    t = cb - R @ ca
    return [R[0, 0], R[0, 1], R[1, 0], R[1, 1], t[0], t[1]], np.hypot(*((R @ A.T).T + t - B).T)


def compose(Mi, M):
    a, b, c, d, e, f = M
    a2, b2, c2, d2, e2, f2 = Mi
    return [a2 * a + b2 * c, a2 * b + b2 * d, c2 * a + d2 * c, c2 * b + d2 * d, a2 * e + b2 * f + e2, c2 * e + d2 * f + f2]


def rot_of(M):
    return math.degrees(math.atan2(M[2], M[0]))


def with_rotation(M, U, deg):
    c = affine_transform(U, M).centroid
    d = math.radians(deg - rot_of(M))
    ct, st = math.cos(d), math.sin(d)
    return compose([ct, -st, st, ct, c.x - (ct * c.x - st * c.y), c.y - (st * c.x + ct * c.y)], M)


def close_onto(gdf, nb_polys, clamp=3.0, maxmove=25.0, Ds=(25, 18, 12, 8, 5, 4, 4, 4)):
    """gdf: the sheet as placed now (EPSG:32644); nb_polys: polygons of the placed neighbours it must close onto"""
    Us = make_valid(unary_union(list(gdf.geometry)))
    B, NB = [], []
    for pp in nb_polys:
        b, nb = dens_pts(pp, 1.0)
        if len(b):
            B.append(b); NB.append(nb)
    if not B:
        return gdf, {"note": "no neighbours"}
    B = np.vstack(B); NB = np.vstack(NB); tB = cKDTree(B)
    M = [1, 0, 0, 1, 0, 0]
    th0 = 0.0
    hist = []

    def matches(M, D):
        P, NP = dens_pts(affine_transform(Us, M), 1.5)
        K = min(40, len(B))
        d, j = tB.query(P, k=K)
        facing = np.einsum("ij,ikj->ik", NP, NB[j]) < -0.3
        ok = (d < D) & facing
        first = np.argmax(ok, axis=1)
        sel = ok[np.arange(len(P)), first]
        return P, j[np.arange(len(P)), first], sel

    for D in Ds:
        P, j, sel = matches(M, D)
        if sel.sum() < 8:
            hist.append((D, int(sel.sum()), None)); continue
        Mi, res = kabsch(P[sel], B[j[sel]])
        M = compose(Mi, M)
        dth = (rot_of(M) - th0 + 180) % 360 - 180
        if abs(dth) > clamp:
            M = with_rotation(M, Us, th0 + math.copysign(clamp, dth))
            P, j, sel = matches(M, D)
            if sel.sum() >= 8:
                t = (B[j[sel]] - P[sel]).mean(0)
                M = compose([1, 0, 0, 1, float(t[0]), float(t[1])], M)
                res = np.hypot(*(P[sel] + t - B[j[sel]]).T)
        hist.append((D, int(sel.sum()), round(float(np.sqrt((res ** 2).mean())), 2)))
    c0, c1 = Us.centroid, affine_transform(Us, M).centroid
    move = math.hypot(c1.x - c0.x, c1.y - c0.y)
    if move > maxmove:
        return gdf, {"note": f"rejected: moved {move:.1f} m > {maxmove}", "hist": hist}
    out = gdf.copy()
    out["geometry"] = [affine_transform(x, M) for x in gdf.geometry]
    return out, {"move_m": round(move, 1), "turn_deg": round(rot_of(M), 2), "hist": hist}
