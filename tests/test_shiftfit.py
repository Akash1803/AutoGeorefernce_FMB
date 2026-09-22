import numpy as np

from autogeoref import shiftfit


def test_survey_keys_from_the_three_sources_normalise_to_one():
    assert shiftfit.normalise("10 B") == shiftfit.normalise("10B") == "10B"
    assert shiftfit.normalise(" 569c ") == "569C"


def test_a_few_control_points_give_one_mean_shift_to_the_whole_village():
    points = np.array([[0.0, 0.0], [100.0, 0.0]])
    disp = np.array([[2.0, -6.0], [4.0, -10.0]])
    shifts, method = shiftfit.fit(points, disp, np.array([[50.0, 50.0], [900.0, 900.0]]))
    assert method == "mean shift"
    # both targets move by the same amount, whatever their distance from the control
    assert np.allclose(shifts[0], [3.0, -8.0]) and np.allclose(shifts[1], [3.0, -8.0])


def test_enough_control_lets_each_parcel_follow_its_neighbours():
    # a village where the error is +10 m north at one end and -10 m at the other
    points = np.array([[0.0, 0.0], [0.0, 100.0], [100.0, 0.0], [100.0, 100.0], [50.0, 50.0]])
    disp = np.array([[0.0, 10.0], [0.0, 10.0], [0.0, -10.0], [0.0, -10.0], [0.0, 0.0]])
    shifts, method = shiftfit.fit(points, disp, np.array([[5.0, 50.0], [95.0, 50.0]]))
    assert method == "local field"
    assert shifts[0][1] > 3.0 and shifts[1][1] < -3.0, "each end keeps its own sign"


def test_a_target_on_top_of_a_control_point_takes_its_displacement_exactly():
    points = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0], [10.0, 10.0]])
    disp = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0]])
    shifts, _ = shiftfit.fit(points, disp, points[2][None, :])
    assert np.allclose(shifts[0], [3.0, 3.0])


def test_leave_one_out_reports_the_error_on_unseen_control():
    # a pure translation: the fit predicts a held-out point perfectly
    points = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0], [100.0, 100.0], [50.0, 50.0]])
    disp = np.repeat(np.array([[3.0, -7.0]]), 5, axis=0)
    loo = shiftfit.leave_one_out(points, disp)
    assert len(loo) == 5 and loo.max() < 1e-6


def test_accuracy_states_both_the_error_before_and_after():
    points = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    disp = np.array([[0.0, -10.0], [0.0, -10.0], [0.0, -10.0]])
    acc = shiftfit.accuracy(points, disp)
    assert acc["control"] == 3
    assert acc["before_median_m"] == 10.0
    assert acc["after_median_m"] == 0.0


def test_no_control_moves_nothing():
    shifts, method = shiftfit.fit(np.zeros((0, 2)), np.zeros((0, 2)), np.array([[1.0, 2.0]]))
    assert method == "none" and np.allclose(shifts[0], [0.0, 0.0])


def test_a_control_parcel_far_from_its_puvi_twin_is_not_control():
    # the same survey number exists in two villages; the wrong one sits kilometres away
    from autogeoref import shift_puvi
    assert shift_puvi.MAX_CONTROL_M <= 200.0, "a plausible cadastral error, not a village apart"
