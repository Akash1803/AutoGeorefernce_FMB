import numpy as np
import pytest

from autogeoref import fit
from autogeoref.match import LineObs, PairObs


def test_a_free_parcel_snaps_onto_a_fixed_anchor_boundary():
    # the free parcel's local edge x = 0 must land on the anchor's ground line x = 100
    line = ((100.0, 0.0), (100.0, 50.0))
    obs = [LineObs("B", (0.0, y), line, 0.30, 0.0) for y in (5.0, 15.0, 25.0, 35.0)]
    out = fit.block_adjust(free={"B": (0.0, np.array([98.0, 0.0]))}, fixed={},
                           pair_obs=[], line_obs=obs, gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (98.0, 0.0), 10.0, 10.0)])
    assert out["B"]["t"][0] == pytest.approx(100.0, abs=0.05)
    assert abs(out["B"]["dtheta_deg"]) < 0.5


def test_a_fixed_anchor_never_moves():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"B": (0.0, np.array([12.0, 1.0]))},
                           fixed={"A": (0.0, np.array([0.0, 0.0]))},
                           pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (12.0, 1.0), 10.0, 10.0)])
    assert "A" not in out
    assert out["B"]["t"] == pytest.approx(np.array([10.0, 0.0]), abs=0.05)


def test_two_free_parcels_meet_each_other():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"A": (0.0, np.array([0.0, 0.0])), "B": (0.0, np.array([11.0, 0.6]))},
                           fixed={}, pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("A", 0.0, (0.0, 0.0), 2.5, 3.0),
                                   fit.PosePrior("B", 0.0, (11.0, 0.6), 2.5, 3.0)])
    a = fit.transform_points([(10.0, 0.0)], out["A"]["theta"], out["A"]["t"])[0]
    b = fit.transform_points([(0.0, 0.0)], out["B"]["theta"], out["B"]["t"])[0]
    assert float(np.linalg.norm(a - b)) < 0.15


def test_gcp_observations_pull_a_parcel_with_no_neighbour():
    gcps = [fit.GcpObs("C", (0.0, 0.0), (500.0, 500.0), 1.0),
            fit.GcpObs("C", (30.0, 0.0), (530.0, 500.0), 1.0)]
    out = fit.block_adjust(free={"C": (0.0, np.array([495.0, 498.0]))}, fixed={},
                           pair_obs=[], line_obs=[], gcp_obs=gcps,
                           priors=[fit.PosePrior("C", 0.0, (495.0, 498.0), 10.0, 10.0)])
    assert out["C"]["t"] == pytest.approx(np.array([500.0, 500.0]), abs=0.2)


def test_the_solver_works_in_raw_utm_coordinates():
    """2026-09-17 bug: with unknowns near 1.4e6 the solver stopped at x0 and moved nothing."""
    line = ((392_500.0, 1_413_000.0), (392_500.0, 1_413_050.0))
    obs = [LineObs("B", (0.0, y), line, 0.30, 0.0) for y in (5.0, 15.0, 25.0, 35.0)]
    start = np.array([392_498.0, 1_413_000.0])
    out = fit.block_adjust(free={"B": (0.0, start)}, fixed={}, pair_obs=[], line_obs=obs, gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, tuple(start), 10.0, 10.0)])
    assert out["B"]["t"][0] == pytest.approx(392_500.0, abs=0.05)
    assert out["B"]["shift_m"] > 1.0


def test_residuals_by_pair_reports_rms_and_max():
    obs = [PairObs("A", "B", (10.0, 0.0), (0.0, 0.0), 0.30),
           PairObs("A", "B", (10.0, 20.0), (0.0, 20.0), 0.30)]
    out = fit.block_adjust(free={"B": (0.0, np.array([12.0, 1.0]))},
                           fixed={"A": (0.0, np.array([0.0, 0.0]))},
                           pair_obs=obs, line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("B", 0.0, (12.0, 1.0), 10.0, 10.0)])
    poses = dict(out)
    poses["A"] = {"theta": 0.0, "t": np.array([0.0, 0.0])}
    r = fit.residuals_by_pair(poses, obs)
    assert r[("A", "B")]["n"] == 2 and r[("A", "B")]["rms"] < 0.15


def test_a_prior_holds_a_parcel_that_has_no_observations():
    out = fit.block_adjust(free={"D": (12.0, np.array([700.0, 800.0]))}, fixed={},
                           pair_obs=[], line_obs=[], gcp_obs=[],
                           priors=[fit.PosePrior("D", 12.0, (700.0, 800.0), 2.5, 3.0)])
    assert out["D"]["shift_m"] < 0.01 and abs(out["D"]["dtheta_deg"]) < 0.01
