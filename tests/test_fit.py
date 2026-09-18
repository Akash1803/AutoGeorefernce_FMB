import numpy as np
import pytest
from shapely import affinity
from shapely.geometry import Polygon

from autogeoref import fit

SQUARE = Polygon([(0, 0), (10, 0), (10, 6), (0, 6)])


def test_rigid_fit_recovers_a_known_rotation_and_shift():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = fit.transform_points(P, 33.0, np.array([150.0, -70.0]))
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    assert theta == pytest.approx(33.0, abs=1e-6)
    assert t == pytest.approx(np.array([150.0, -70.0]), abs=1e-6)
    assert rms < 1e-9 and mx < 1e-9


def test_rigid_fit_never_scales_even_when_the_target_is_stretched():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P * 1.10
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    moved = fit.transform_points(P, theta, t)
    side = float(np.linalg.norm(moved[1] - moved[0]))
    assert side == pytest.approx(10.0, abs=1e-9), "scale must stay 1"
    assert rms > 0.1, "the stretch shows up as residual, not as scale"


def test_rigid_fit_rejects_a_mirror():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P.copy()
    Q[:, 1] *= -1
    theta, t, rms, mx = fit.rigid_fit(P, Q)
    moved = fit.transform_points(P, theta, t)
    u, v = moved[1] - moved[0], moved[2] - moved[1]
    assert u[0] * v[1] - u[1] * v[0] > 0, "handedness preserved"   # 2D cross product


def test_similarity_fit_reports_scale_and_anisotropy():
    P = np.array([(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)])
    Q = P * 0.90
    theta, t, rms, scale, aniso = fit.similarity_fit(P, Q)
    assert scale == pytest.approx(0.90, abs=1e-6)
    assert aniso == pytest.approx(1.0, abs=1e-6)


def test_apply_pose_preserves_area_and_perimeter():
    moved = fit.apply_pose(SQUARE, 47.5, np.array([1000.0, 2000.0]))
    assert moved.area == pytest.approx(SQUARE.area, abs=1e-9)
    assert moved.length == pytest.approx(SQUARE.length, abs=1e-9)
    fit.assert_rigid(SQUARE, moved)


def test_assert_rigid_catches_a_stretch():
    with pytest.raises(AssertionError):
        fit.assert_rigid(SQUARE, affinity.scale(SQUARE, 1.01, 1.01))
