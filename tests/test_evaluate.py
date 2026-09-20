import pytest

from autogeoref import evaluate, paths

VILLAGE = "35_04_077"


def test_a_loose_points_file_falls_back_to_the_hand_geometry():
    """47B's three control points fit at 6.4 m rms; its placed geometry fits at 1.8 m."""
    if not paths.vector_dir(VILLAGE).exists():
        pytest.skip("village not present")
    from autogeoref import anchors
    assert anchors.pose_from_points(VILLAGE, "47B")[2] > evaluate.REPORT_ONLY_RMS
    t = evaluate.truth_pose(VILLAGE, "47B")
    assert t["source"] == "geometry"
    assert t["rms"] < evaluate.REPORT_ONLY_RMS


def test_the_heading_spread_needs_four_points_to_mean_anything():
    """Drop one of three and the remaining two fit any heading exactly, so the spread is 0."""
    if not paths.vector_dir(VILLAGE).exists():
        pytest.skip("village not present")
    assert evaluate.truth_pose(VILLAGE, "47B")["heading_spread"] == 0.0     # 3 points
    assert evaluate.truth_pose(VILLAGE, "48A")["heading_spread"] > 0.0      # 4 points


def test_tolerances_never_tighter_than_the_truth_itself():
    m, deg = evaluate.tolerances({"rms": 5.85, "heading_spread": 1.6})
    assert m == pytest.approx(5.85) and deg == pytest.approx(1.6)
    m2, deg2 = evaluate.tolerances({"rms": 0.19, "heading_spread": 0.3})
    assert m2 == pytest.approx(1.5) and deg2 == pytest.approx(1.5)


def test_calibrate_separates_truth_from_the_best_wrong_pose():
    rows = [{"survey": "48A", "share_truth": 0.72, "share_best_wrong": 0.44},
            {"survey": "47B", "share_truth": 0.48, "share_best_wrong": 0.30},
            {"survey": "46B", "share_truth": 0.12, "share_best_wrong": 0.46}]
    cal = evaluate.calibrate(rows)
    assert 0.0 < cal["amber_share"] < cal["green_share"] <= 1.0
    assert cal["green_share"] > 0.46, "a threshold below the best wrong pose would pass 46B wrongly"


def test_the_harness_measures_real_ground_error_not_the_solver_shift():
    """shift_m is movement from the start guess; the score must be distance from the team's truth."""
    import numpy as np
    from shapely.geometry import Polygon
    truth = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    got = Polygon([(3, 4), (13, 4), (13, 14), (3, 14)])
    assert evaluate.centroid_error(got, truth) == pytest.approx(5.0, abs=1e-9)
    assert evaluate.heading_error(359.0, 1.0) == pytest.approx(2.0, abs=1e-9)
    assert np.isfinite(evaluate.centroid_error(got, truth))


@pytest.mark.slow
def test_leave_one_out_runs_on_a_copy_and_never_touches_the_project():
    if not paths.vector_dir(VILLAGE).exists():
        pytest.skip("village not present")
    before = {p.name: p.stat().st_mtime
              for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    out = evaluate.leave_one_out(VILLAGE, surveys=["43A"], mode="full")
    after = {p.name: p.stat().st_mtime
             for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    assert before == after, "the harness must work on a copy"
    assert out and set(out[0]) >= {"survey", "colour", "centroid_error_m", "heading_error_deg",
                                   "passed"}


def test_ring_around_a_seed_uses_the_printed_neighbours_both_ways():
    if not paths.vector_dir(VILLAGE).exists():
        pytest.skip("village not present")
    ring1 = evaluate.ring_around(VILLAGE, ["48A"], rings=1)
    assert "47B" in ring1 and "171" in ring1 and "48A" not in ring1
    ring2 = evaluate.ring_around(VILLAGE, ["48A"], rings=2)
    assert set(ring1) < set(ring2)
