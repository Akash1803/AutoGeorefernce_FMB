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
