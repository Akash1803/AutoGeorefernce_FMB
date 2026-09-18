"""Rigid (rotation + shift, scale 1) fitting and pose application."""
import math

import numpy as np
from shapely import affinity


def _procrustes(P, Q):
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:                      # reflection: flip the smallest singular direction
        Vt = Vt.copy()
        Vt[-1] *= -1
        R = Vt.T @ U.T
    return R, pc, qc, S, P, Q


def rigid_fit(P, Q):
    """Rotation + shift with scale fixed at 1. Returns (theta_deg, t, rms, max_residual)."""
    R, pc, qc, _S, P, Q = _procrustes(P, Q)
    t = qc - R @ pc
    res = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
    theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    return theta, t, float(np.sqrt((res ** 2).mean())), float(res.max())


def similarity_fit(P, Q):
    """Descriptive only: how much scale and anisotropy a hand placement carries."""
    R, pc, qc, S, P, Q = _procrustes(P, Q)
    denom = ((P - pc) ** 2).sum()
    scale = float(S.sum() / denom) if denom else float("nan")
    t = qc - scale * (R @ pc)
    res = np.linalg.norm(scale * (R @ P.T).T + t - Q, axis=1)
    A = np.hstack([P, np.ones((len(P), 1))])
    M, *_ = np.linalg.lstsq(A, Q, rcond=None)
    sv = np.linalg.svd(M[:2].T, compute_uv=False)
    aniso = float(sv.max() / sv.min()) if sv.min() > 0 else float("inf")
    theta = math.degrees(math.atan2(R[1, 0], R[0, 0]))
    return theta, t, float(np.sqrt((res ** 2).mean())), scale, aniso


def transform_points(P, theta_deg, t):
    th = math.radians(theta_deg)
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    return (R @ np.asarray(P, float).T).T + np.asarray(t, float)


def apply_pose(geom, theta_deg, t):
    return affinity.translate(affinity.rotate(geom, theta_deg, origin=(0, 0)), float(t[0]), float(t[1]))


def assert_rigid(before, after, tol=1e-6):
    if abs(after.area - before.area) > tol or abs(after.length - before.length) > tol:
        raise AssertionError(
            "placement is not rigid: area %.9f -> %.9f, perimeter %.9f -> %.9f"
            % (before.area, after.area, before.length, after.length))


# --------------------------------------------------------------------------- block adjustment
from dataclasses import dataclass as _dataclass                 # noqa: E402

from scipy.optimize import least_squares                        # noqa: E402


@_dataclass
class GcpObs:
    """A sheet point that should land on a known ground point."""
    survey: str
    p_local: tuple
    map_xy: tuple
    sigma: float


@_dataclass
class PosePrior:
    """Where a parcel is believed to be before the adjustment, and how strongly."""
    survey: str
    theta: float
    t: tuple
    sigma_pos: float
    sigma_head: float


def block_adjust(free, fixed, pair_obs, line_obs, gcp_obs, priors, f_scale=1.5):
    """Least-squares rotation + shift per free survey; `fixed` surveys never move.

    Solved in a frame centred on the block: with raw UTM unknowns (easting about 3.9e5, northing
    about 1.4e6) scipy's least_squares terminates on xtol at the starting point and nothing moves.
    That was the 2026-09-17 bug.
    """
    names = sorted(free)
    if not names:
        return {}
    idx = {s: i for i, s in enumerate(names)}
    origin = np.mean([np.asarray(t, float) for _th, t in free.values()], axis=0)

    x0 = np.zeros(3 * len(names))
    for s in names:
        th, t = free[s]
        x0[3 * idx[s]] = math.radians(th)
        x0[3 * idx[s] + 1:3 * idx[s] + 3] = np.asarray(t, float) - origin

    fixed_pose = {s: (math.radians(th), np.asarray(t, float) - origin) for s, (th, t) in fixed.items()}

    def place(survey, x, p):
        p = np.asarray(p, float)
        if survey in fixed_pose:
            th, t = fixed_pose[survey]
        else:
            i = idx[survey]
            th, t = x[3 * i], x[3 * i + 1:3 * i + 3]
        c, s = math.cos(th), math.sin(th)
        return np.array([c * p[0] - s * p[1] + t[0], s * p[0] + c * p[1] + t[1]])

    def residuals(x):
        r = []
        for o in pair_obs:
            r += list((place(o.a, x, o.pa) - place(o.b, x, o.pb)) / o.sigma)
        for o in line_obs:
            p = place(o.survey, x, o.p_local)
            a = np.asarray(o.line[0], float) - origin
            b = np.asarray(o.line[1], float) - origin
            d = b - a
            n = float(np.hypot(d[0], d[1]))
            perp = ((p[0] - a[0]) * d[1] - (p[1] - a[1]) * d[0]) / n if n else 0.0
            r.append(perp / o.sigma)
        for o in gcp_obs:
            r += list((place(o.survey, x, o.p_local) - (np.asarray(o.map_xy, float) - origin)) / o.sigma)
        for pr in priors:
            if pr.survey not in idx:
                continue
            i = idx[pr.survey]
            r += list((x[3 * i + 1:3 * i + 3] - (np.asarray(pr.t, float) - origin)) / pr.sigma_pos)
            d = (x[3 * i] - math.radians(pr.theta) + math.pi) % (2 * math.pi) - math.pi
            r.append(d / math.radians(pr.sigma_head))
        return np.array(r) if r else np.zeros(1)

    sol = least_squares(residuals, x0, loss="huber", f_scale=f_scale, x_scale="jac",
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    try:
        cov = np.linalg.pinv(sol.jac.T @ sol.jac)
        sig = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        sig = np.full(3 * len(names), np.nan)

    out = {}
    for s in names:
        i = idx[s]
        th = math.degrees(sol.x[3 * i])
        t = sol.x[3 * i + 1:3 * i + 3] + origin
        th0, t0 = free[s]
        out[s] = {"theta": float(th), "t": t,
                  "dtheta_deg": float(((th - th0) + 180) % 360 - 180),
                  "shift_m": float(np.linalg.norm(t - np.asarray(t0, float))),
                  "sigma_head_deg": float(math.degrees(sig[3 * i])),
                  "sigma_pos_m": float(np.hypot(sig[3 * i + 1], sig[3 * i + 2]))}
    return out


def residuals_by_pair(poses, pair_obs):
    """{(a, b): {'n', 'rms', 'max'}} after an adjustment; `poses` maps survey -> {'theta', 't'}."""
    acc = {}
    for o in pair_obs:
        if o.a not in poses or o.b not in poses:
            continue
        pa = transform_points([o.pa], poses[o.a]["theta"], poses[o.a]["t"])[0]
        pb = transform_points([o.pb], poses[o.b]["theta"], poses[o.b]["t"])[0]
        acc.setdefault((o.a, o.b), []).append(float(np.linalg.norm(pa - pb)))
    return {k: {"n": len(v), "rms": float(np.sqrt(np.mean(np.square(v)))), "max": float(max(v))}
            for k, v in acc.items()}
