import numpy as np

from autogeoref import shift_score


def test_rmse_is_the_root_mean_square_not_the_mean():
    assert abs(shift_score.rmse([3.0, 4.0]) - 3.5355) < 0.001
    assert shift_score.rmse([]) != shift_score.rmse([])          # nan on an empty set


def test_summarise_counts_what_improved_and_what_did_not():
    rows = [{"err_puvi_m": 10.0, "err_layer_m": 2.0},
            {"err_puvi_m": 4.0, "err_layer_m": 9.0},
            {"err_puvi_m": 5.0, "err_layer_m": 5.0}]
    s = shift_score.summarise(rows)
    assert s["parcels"] == 3 and s["better"] == 1 and s["worse"] == 1
    assert s["within_3m"] == 1 and s["within_10m"] == 3
    assert s["layer_median_m"] == 5.0


def test_a_control_parcel_a_village_away_is_not_truth():
    # the guard that stopped Thailavaram parcels being used as control for Potheri
    from autogeoref import shift_puvi
    assert shift_puvi.MAX_CONTROL_M <= 200.0
