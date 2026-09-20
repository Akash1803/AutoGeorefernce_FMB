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


def test_a_parcel_is_turned_onto_a_strips_edge_by_line_observations():
    """47B / 171, 2026-09-20: the strip draws one 300 m edge, the parcel its 80 m share of it,
    so no equal-length chain exists; points of the parcel's edge must still lie on the line."""
    import math
    strip_local = [(0.0, 0.0), (300.0, 0.0), (300.0, 15.0), (0.0, 15.0)]      # fixed, at identity
    parcel_local = [(0.0, 0.0), (80.0, 0.0), (80.0, 40.0), (0.0, 40.0)]
    # the parcel belongs at (100, 15) on top of the strip; start it 6 degrees off about its corner
    truth_theta, truth_t = 0.0, np.array([100.0, 15.0])
    start_theta = 6.0
    obs = []
    for k in range(0, 81, 2):                                                   # samples along its bottom edge
        obs.append(fit.PairLineObs("P", (float(k), 0.0), "S", (0.0, 15.0), (300.0, 15.0), 0.6))
    # one short chain pins the corner (like 47B's short edge with 48A)
    from autogeoref import match
    pair = [match.PairObs("P", "S", (0.0, 0.0), (100.0, 15.0), 0.3)]
    prior = [fit.PosePrior("P", start_theta, tuple(truth_t), 10.0, 10.0)]
    out = fit.block_adjust({"P": (start_theta, truth_t.copy())}, {"S": (0.0, np.zeros(2))},
                           pair, [], [], prior, pair_line_obs=obs)
    assert abs(out["P"]["theta"] - truth_theta) < 0.3, out["P"]["theta"]
