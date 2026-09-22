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
    # a village a kilometre across where the error is +10 m north at one end, -10 m at the other
    points = np.array([[0.0, 0.0], [0.0, 1000.0], [1000.0, 0.0], [1000.0, 1000.0], [500.0, 500.0]])
    disp = np.array([[0.0, 10.0], [0.0, 10.0], [0.0, -10.0], [0.0, -10.0], [0.0, 0.0]])
    shifts, method = shiftfit.fit(points, disp, np.array([[50.0, 500.0], [950.0, 500.0]]))
    assert method == "local field"
    # the field is deliberately gentle, so the ends are pulled apart but not all the way to the
    # +-10 m the control says: what matters is that each end follows the control nearest it
    assert shifts[0][1] > 2.0 and shifts[1][1] < -2.0, "each end keeps its own sign"
    assert abs(shifts[0][1] + shifts[1][1]) < 0.5, "and symmetrically"


def test_a_target_leans_towards_the_control_it_sits_on_without_being_yanked_onto_it():
    # the field is smoothed on purpose: a control point a metre away must not outvote the rest,
    # or a parcel straddling it is torn. So the answer leans to 3.0 without reaching it.
    points = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0], [10.0, 10.0]])
    disp = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0]])
    shifts, _ = shiftfit.fit(points, disp, points[2][None, :])
    others = shiftfit.fit(points, disp, np.array([[5.0, 5.0]]))[0][0]
    assert shifts[0][0] > others[0], "closer to the control it sits on than the middle of the village is"
    assert 2.0 < shifts[0][0] < 3.0


def test_the_smoothing_distance_is_what_keeps_a_parcel_from_being_torn():
    # two vertices of one parcel, a metre apart, must not be given wildly different moves
    points = np.array([[0.0, 0.0], [300.0, 0.0], [0.0, 300.0], [300.0, 300.0], [150.0, 150.0]])
    disp = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [30.0, 0.0]])
    a = shiftfit.idw(points, disp, np.array([150.0, 150.0]))
    b = shiftfit.idw(points, disp, np.array([151.0, 150.0]))
    assert abs(a[0] - b[0]) < 0.2, "a metre apart, moved within 20 cm of each other"


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


def test_two_parcels_that_share_a_boundary_still_share_it_after_the_warp():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    # the field pulls harder on the right than on the left
    points = np.array([[0.0, 0.0], [0.0, 400.0], [800.0, 0.0], [800.0, 400.0]])
    disp = np.array([[0.0, 0.0], [0.0, 0.0], [20.0, 5.0], [20.0, 5.0]])
    left = Polygon([(0, 0), (400, 0), (400, 400), (0, 400)])
    right = Polygon([(400, 0), (800, 0), (800, 400), (400, 400)])
    table = shift_puvi._warp_table(shift_puvi._nodes([left, right]), points, disp)
    a = shift_puvi._warp_with(left, table)
    b = shift_puvi._warp_with(right, table)
    assert a.intersection(b).area < 0.01, "no sliver of overlap"
    merged = a.union(b)
    assert merged.geom_type == "Polygon" and len(list(merged.interiors)) == 0, "and no gap between them"
    assert a.boundary.buffer(0.01).intersection(b.boundary).length > 380, "the shared line is one line"


def test_the_warp_moves_a_parcel_by_the_field():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    points = np.array([[0.0, 0.0], [500.0, 0.0], [0.0, 500.0], [500.0, 500.0]])
    disp = np.repeat(np.array([[3.0, -7.0]]), 4, axis=0)
    p = Polygon([(100, 100), (300, 100), (300, 300), (100, 300)])
    w = shift_puvi._warp_with(p, shift_puvi._warp_table(shift_puvi._nodes([p]), points, disp))
    assert abs(w.centroid.x - (p.centroid.x + 3.0)) < 0.01
    assert abs(w.centroid.y - (p.centroid.y - 7.0)) < 0.01
    assert abs(w.area - p.area) < 0.5, "a uniform field is a plain shift: the area is unchanged"


def test_neighbours_whose_vertices_differ_by_millimetres_still_move_as_one():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    points = np.array([[0.0, 0.0], [0.0, 400.0], [800.0, 0.0], [800.0, 400.0]])
    disp = np.array([[0.0, 0.0], [0.0, 0.0], [20.0, 5.0], [20.0, 5.0]])
    left = Polygon([(0, 0), (400, 0), (400, 400), (0, 400)])
    right = Polygon([(400.003, 0.002), (800, 0), (800, 400), (399.998, 399.997)])  # as Puvi draws it
    table = shift_puvi._warp_table(shift_puvi._nodes([left, right]), points, disp)
    a, b = shift_puvi._warp_with(left, table), shift_puvi._warp_with(right, table)
    assert a.intersection(b).area < 0.001, "no hairline overlap"
    merged = a.union(b)
    assert merged.geom_type == "Polygon" and len(list(merged.interiors)) == 0, "no hairline gap"


def test_seam_control_pulls_a_village_edge_onto_its_neighbour():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    # two villages digitised apart: a 20 m gap along a 400 m frontage
    left = Polygon([(0, 0), (400, 0), (400, 400), (0, 400)])
    right = Polygon([(420, 0), (800, 0), (800, 400), (420, 400)])
    pts, disp = shift_puvi.seam_control(right, left, share=1.0)
    assert len(pts) > 5, "the frontage is sampled"
    facing = [d for p_, d in zip(pts, disp) if p_[0] < 430]
    assert facing and all(d[0] < -5.0 for d in facing), "the facing edge is pulled west onto the neighbour"


def test_seam_control_shares_the_move_when_neither_side_is_an_authority():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    left = Polygon([(0, 0), (400, 0), (400, 400), (0, 400)])
    right = Polygon([(420, 0), (800, 0), (800, 400), (420, 400)])
    _, full = shift_puvi.seam_control(right, left, share=1.0)
    _, half = shift_puvi.seam_control(right, left, share=0.5)
    assert abs(half[0][0] * 2 - full[0][0]) < 1e-6, "half the distance each"


def test_seam_control_is_quiet_when_the_villages_already_meet():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    left = Polygon([(0, 0), (400, 0), (400, 400), (0, 400)])
    right = Polygon([(400, 0), (800, 0), (800, 400), (400, 400)])
    pts, disp = shift_puvi.seam_control(right, left)
    moves = [abs(d[0]) + abs(d[1]) for p_, d in zip(pts, disp) if p_[0] < 410]
    assert not moves or max(moves) < 1.0, "nothing to close along a boundary already shared"


def test_a_village_seam_overlap_is_clipped_out_of_the_side_without_our_placements():
    import geopandas as gpd
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    anchored = gpd.GeoDataFrame(geometry=[Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])], crs=32644)
    loose = gpd.GeoDataFrame(geometry=[Polygon([(90, 0), (200, 0), (200, 100), (90, 100)])], crs=32644)
    layers, notes = shift_puvi.clip_village_overlaps({"A": anchored, "B": loose}, {"A": True, "B": False})
    assert layers["A"].geometry.iloc[0].area == 10000, "the anchored village is untouched"
    assert abs(layers["B"].geometry.iloc[0].area - 10000) < 1.0, "the other gives up the disputed strip"
    assert notes and "B yielded" in notes[0]


def test_the_clip_refuses_to_eat_a_parcel():
    import geopandas as gpd
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    anchored = gpd.GeoDataFrame(geometry=[Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])], crs=32644)
    # this parcel lies almost entirely inside the anchored village: moving it is a placement
    # question, not a seam, so the guard leaves it alone
    loose = gpd.GeoDataFrame(geometry=[Polygon([(10, 10), (90, 10), (90, 90), (10, 90)])], crs=32644)
    layers, notes = shift_puvi.clip_village_overlaps({"A": anchored, "B": loose}, {"A": True, "B": False})
    assert layers["B"].geometry.iloc[0].area == 6400, "untouched"
    assert "refused by the guard" in notes[0]


def test_the_reach_beyond_which_a_correction_is_not_evidenced_is_the_measured_one():
    from autogeoref import shift_puvi
    # measured on 2026-09-22: the fit worked under 150 m and did nothing at 150-400 m
    assert shift_puvi.MEASURED_REACH_M == 150.0


def test_two_parcels_of_one_village_that_overlap_are_settled_by_the_larger_one():
    import geopandas as gpd
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    big = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
    small = Polygon([(95, 0), (150, 0), (150, 100), (95, 100)])
    gdf = gpd.GeoDataFrame(geometry=[big, small], crs=32644)
    out, fixed = shift_puvi.clip_siblings(gdf)
    assert fixed == 1
    assert out.geometry.iloc[0].intersection(out.geometry.iloc[1]).area < 0.01
    assert abs(out.geometry.iloc[1].area - 5500) < 1.0, "the smaller parcel keeps its ground"


def test_a_parcel_stored_as_a_mixed_collection_is_still_warped():
    from shapely.geometry import GeometryCollection, LineString, Polygon
    from autogeoref import shift_puvi
    # Puvi stores Peramanur survey 11 like this: a polygon and a stray line in one geometry
    poly = Polygon([(10, 10), (60, 10), (60, 60), (10, 60)])
    mixed = GeometryCollection([poly, LineString([(0, 0), (5, 5)])])
    points = np.array([[0.0, 0.0], [500.0, 0.0], [0.0, 500.0], [500.0, 500.0]])
    disp = np.repeat(np.array([[4.0, -3.0]]), 4, axis=0)
    table = shift_puvi._warp_table(shift_puvi._nodes([mixed]), points, disp)
    out = shift_puvi._warp_with(mixed, table)
    assert out is not None and not out.is_empty, "the parcel survives"
    assert abs(out.area - poly.area) < 1.0
    assert abs(out.centroid.x - (poly.centroid.x + 4.0)) < 0.2


def test_a_parcel_is_never_lost_even_if_the_warp_fails():
    from shapely.geometry import Polygon
    from autogeoref import shift_puvi
    poly = Polygon([(10, 10), (60, 10), (60, 60), (10, 60)])
    out = shift_puvi._warp_with(poly, {}, field=lambda q: np.array([5.0, 0.0]))
    assert out is not None and not out.is_empty


def test_the_corridor_list_is_every_buffer_village_that_has_puvi_data():
    from autogeoref import shift_puvi
    # 23 revenue villages are crossed by the 30 m buffer; Tambaram 35_05_010 has no Puvi vector
    assert len(shift_puvi.CORRIDOR) == 22
    assert "35_05_010" not in shift_puvi.CORRIDOR
    assert set(shift_puvi.STRETCH) <= set(shift_puvi.CORRIDOR)
    assert len(set(shift_puvi.CORRIDOR)) == len(shift_puvi.CORRIDOR), "no village listed twice"


def test_a_seam_point_is_not_evidence():
    import inspect
    from autogeoref import shift_puvi
    src = inspect.getsource(shift_puvi.run_village)
    assert "hand_points = cpoints.copy()" in src
    assert "np.linalg.norm(hand_points - q" in src, "the evidence flag measures distance to the team's own parcels"


def test_a_crowd_of_seam_points_is_averaged_more_widely_than_hand_control():
    # 40 seam points along one frontage must not yank the parcels beside them out of shape
    seam = np.array([[float(i) * 10.0, 0.0] for i in range(40)])
    disp = np.repeat(np.array([[0.0, 8.0]]), 40, axis=0)
    disp[20:] = [0.0, -8.0]                       # the frontage disagrees with itself mid-way
    q = np.array([[195.0, 30.0], [205.0, 30.0]])
    sharp, _ = shiftfit.fit(seam, disp, q)
    wide, _ = shiftfit.fit(seam, disp, q, k=shiftfit.SEAM_K, smooth=shiftfit.SEAM_SMOOTH_M)
    assert abs(wide[0][1] - wide[1][1]) < abs(sharp[0][1] - sharp[1][1]), "the wide field varies less"
    assert abs(wide[0][1] - wide[1][1]) < 1.0, "two points 10 m apart move together"


def test_one_control_parcel_that_disagrees_with_its_neighbours_stops_dragging_the_answer():
    # eight neighbours agree the ground moves 10 m north; one says 60 m south
    points = np.array([[float(i) * 40.0, 0.0] for i in range(9)])
    disp = np.repeat(np.array([[0.0, 10.0]]), 9, axis=0)
    disp[4] = [0.0, -60.0]
    q = np.array([160.0, 5.0])
    plain = shiftfit.idw(points, disp, q, rounds=0)
    robust = shiftfit.idw(points, disp, q)
    assert abs(robust[1] - 10.0) < abs(plain[1] - 10.0), "the odd one out counts for less"
    assert robust[1] > 0, "and no longer flips the direction"


def test_the_robust_step_leaves_agreeing_control_alone():
    points = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0], [100.0, 100.0], [50.0, 50.0]])
    disp = np.repeat(np.array([[3.0, -7.0]]), 5, axis=0)
    assert np.allclose(shiftfit.idw(points, disp, np.array([50.0, 20.0])), [3.0, -7.0], atol=1e-6)
