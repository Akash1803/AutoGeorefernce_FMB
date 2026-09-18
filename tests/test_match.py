import numpy as np
import pytest
from shapely.geometry import LineString, Polygon

from autogeoref import fit, match, sheets

A = Polygon([(0, 0), (40, 0), (40, 25), (0, 25)])       # shares the x = 40 edge with B
B = Polygon([(40, 0), (75, 0), (75, 25), (40, 25)])


def test_common_chains_finds_the_shared_edge_of_two_rectangles():
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, B)]))
    assert chains, "the 25 m shared edge must be found"
    assert chains[0].length == pytest.approx(25.0, abs=0.5)
    assert chains[0].n_pairs >= 2


def test_common_chains_returns_every_run_not_only_the_longest():
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, B)]))
    assert len(chains) >= 2, "returning only the longest run hid true matches in the prototype"
    assert chains[0].length >= chains[-1].length


def test_a_chain_pins_the_pose_of_a_moved_neighbour():
    moved = fit.apply_pose(B, 20.0, np.array([300.0, 400.0]))
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, moved)]))
    best = chains[0]
    P = np.array([q for _p, q in best.pairs])
    Q = np.array([p for p, _q in best.pairs])
    theta, t, rms, _mx = fit.rigid_fit(P, Q)
    assert rms < 0.05
    assert abs(((theta + 20.0) + 180) % 360 - 180) < 0.5


def test_no_chain_when_no_run_of_edges_has_matching_lengths():
    odd = Polygon([(0, 0), (13.7, 0), (13.7, 9.1), (0, 9.1)])
    assert match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, odd)])) == []


def test_the_matcher_is_position_blind_so_the_caller_must_gate_on_distance():
    """Two congruent parcels match wherever they sit: edge lengths carry no position.

    The engine therefore only offers pairs whose placed geometry is already within a few metres.
    """
    far = Polygon([(500, 500), (540, 500), (540, 525), (500, 525)])   # congruent with A
    assert match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, far)]))


def test_line_observations_measure_perpendicular_distance_to_an_anchor_edge():
    obs = match.line_observations("42B", np.array([[40.5, 5.0], [40.5, 15.0]]),
                                  LineString([(40, 0), (40, 25)]), sigma=0.30)
    assert len(obs) == 2
    assert all(abs(o.distance_now - 0.5) < 1e-9 for o in obs)
    assert all(o.sigma == 0.30 for o in obs)


def test_line_observations_ignore_samples_further_than_max_dist():
    obs = match.line_observations("42B", np.array([[40.5, 5.0], [55.0, 5.0]]),
                                  LineString([(40, 0), (40, 25)]), sigma=0.30, max_dist=2.0)
    assert len(obs) == 1


def test_chain_observations_carry_both_survey_ids_and_sigma():
    chains = match.common_chains(sheets.outline([(None, A)]), sheets.outline([(None, B)]))
    obs = match.chain_observations(chains[0], "A", "B", sigma=0.30)
    assert len(obs) == chains[0].n_pairs
    assert all(o.a == "A" and o.b == "B" and o.sigma == 0.30 for o in obs)
