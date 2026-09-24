import numpy as np
import pytest

from autogeoref import engine


def test_colour_rule_green_needs_an_anchor_and_a_decisive_majority_or_clear_imagery():
    base = {"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.65, "observable": True,
            "n_neighbours": 2, "method": "anchors", "ambiguous": False,
            "support_m": 200.0, "support_next_m": 40.0}
    assert engine.colour_of(base) == "amber", "two neighbours and a middling image: not enough"
    assert engine.colour_of(dict(base, n_neighbours=3)) == "green"
    assert engine.colour_of(dict(base, share=0.90)) == "green"


def test_colour_rule_amber_when_only_one_condition_is_short():
    assert engine.colour_of({"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.20,
                             "observable": False, "n_neighbours": 2, "method": "anchors",
                             "ambiguous": False, "support_m": 100.0,
                             "support_next_m": 0.0}) == "amber"


def test_colour_rule_red_for_an_ambiguous_image_pose():
    assert engine.colour_of({"boundary_rms_m": None, "anchor_rms_m": None, "share": 0.55,
                             "observable": True, "n_neighbours": 0, "method": "image",
                             "ambiguous": True}) == "red"


def test_colour_rule_red_for_puvi_only():
    assert engine.colour_of({"boundary_rms_m": None, "anchor_rms_m": None, "share": 0.0,
                             "observable": False, "n_neighbours": 0, "method": "puvi-only",
                             "ambiguous": True}) == "red"


def test_write_parcels_keeps_fmb_dimensions_exactly(tmp_path):
    import geopandas as gpd
    from autogeoref import paths, sheets
    if not paths.sheet_path("35_04_077", "46B").exists():
        pytest.skip("village sheets not present")
    out = engine.write_parcels("35_04_077", "46B", 31.0, np.array([392800.0, 1413300.0]),
                               {"colour": "green", "method": "anchors", "confidence": 90},
                               out_dir=tmp_path)
    got = gpd.read_file(out, layer="parcels")
    src = sheets.load_sheet("35_04_077", "46B")
    assert len(got) == len(src)
    assert float(got.geometry.area.sum()) == pytest.approx(sum(g.area for _p, g in src), abs=1e-6)
    edges = gpd.read_file(out, layer="edges")
    assert {"length_m", "bearing_grid_deg", "bearing_true_deg"} <= set(edges.columns)
    assert float(edges["length_m"].max()) > 1.0


def test_run_on_a_village_where_every_parcel_is_an_anchor_changes_nothing(village_copy):
    """Kizhikaranai is fully hand placed: the engine must place nothing and touch no manual file."""
    from autogeoref import paths
    village_copy("35_04_077")
    before = {p.name: (p.stat().st_mtime, p.stat().st_size)
              for p in paths.vector_dir("35_04_077").glob("*_parcels_modified*.gpkg")}
    out = engine.run("35_04_077", do_raster=False, do_topology=False, do_review=False)
    after = {p.name: (p.stat().st_mtime, p.stat().st_size)
             for p in paths.vector_dir("35_04_077").glob("*_parcels_modified*.gpkg")}
    assert before == after, "anchors are frozen"
    assert out["anchors"] == 15 and out["placed"] == 0


def test_observations_reach_the_solver_in_sheet_metres():
    """On 2026-09-19 ground points were handed to the solver as if they were sheet points and a
    correctly chosen 46B was carried 72 m. A parcel whose chains already coincide must not move."""
    from autogeoref import fit, paths
    if not paths.vector_dir("35_04_077").exists():
        pytest.skip("village not present")
    from autogeoref import anchors, evaluate
    v, s = "35_04_077", "46B"
    amap = {k: a for k, a in anchors.load_anchors(v, set()).items() if k != s}
    placed = {k: (a.theta, a.t) for k, a in amap.items()}
    truth = evaluate.truth_pose(v, s)
    pose = (truth["theta"], truth["t"])
    obs, _sup, _part = engine.neighbour_chains(v, s, pose, placed, amap)
    assert obs, "46B shares boundary with 43A, 47A and 170"
    for o in obs:
        assert abs(o.pa[0]) < 2000 and abs(o.pa[1]) < 2000, "sheet metres, not UTM"
    adj = fit.block_adjust({s: pose}, placed, obs, [], [],
                           [fit.PosePrior(s, pose[0], tuple(pose[1]), engine.SIGMA_START, 10.0)])
    before = evaluate._sheet_at(v, s, *pose)
    after = evaluate._sheet_at(v, s, adj[s]["theta"], adj[s]["t"])
    assert before.centroid.distance(after.centroid) < 3.0


def test_dedupe_collapses_the_same_pose_from_many_chains():
    base = {"theta": -23.9, "t": np.array([392900.0, 1413500.0]), "from": "170", "chain_m": 70.0,
            "fit_rms": 0.1}
    # the same landing from another chain: rotate by 0.3 degrees about the parcel itself
    from autogeoref import fit
    probe = np.array([[50.0, 50.0]])
    land = fit.transform_points(probe, base["theta"], base["t"])[0]
    th2 = base["theta"] + 0.3
    t2 = land - fit.transform_points(probe, th2, np.zeros(2))[0]
    near = dict(base, theta=th2, t=t2, chain_m=20.0)
    far = dict(base, theta=156.1, chain_m=40.0)                 # the half-turn twin
    out = engine.dedupe_poses([base, near, far])
    assert len(out) == 2
    assert max(c["votes"] for c in out) == 2


def test_a_tie_between_poses_is_never_green():
    """42A, 2026-09-19: 348 m vs 347 m of support, residual 1.5 m, parcel 6.4 m out."""
    tie = {"boundary_rms_m": 1.48, "anchor_rms_m": 2.2, "share": 0.60, "observable": False,
           "n_neighbours": 4, "method": "neighbour", "ambiguous": False,
           "support_m": 348.0, "support_next_m": 347.0}
    assert engine.colour_of(tie) == "amber"
    clear = dict(tie, support_m=271.0, support_next_m=29.0)
    assert engine.colour_of(clear) == "green"
    imaged = dict(tie, share=0.90, observable=True)
    assert engine.colour_of(imaged) == "green", "clear imagery may break a geometric tie"


def test_a_pose_that_misses_a_printed_neighbour_is_red():
    """47A, 2026-09-19: 129 m out, residual 0.75 m, and none of 46B/47B/48A within reach."""
    row = {"boundary_rms_m": 0.75, "anchor_rms_m": 2.2, "share": 0.6, "observable": True,
           "n_neighbours": 1, "method": "neighbour", "ambiguous": False,
           "support_m": 120.0, "support_next_m": 20.0, "printed_far": "46B,48A"}
    assert engine.colour_of(row) == "red"
    assert engine.colour_of(dict(row, printed_far="", n_neighbours=3)) == "green"
    # 43B: its sheet names 42B, but the team placed 42B and 43B 20 m apart. One miss is a doubt.
    assert engine.colour_of(dict(row, printed_far="42B", n_neighbours=3)) == "amber"


def test_printed_contradictions_reads_the_transcription():
    from autogeoref import anchors, evaluate, paths
    if not paths.vector_dir("35_04_077").exists():
        pytest.skip("village not present")
    v, s = "35_04_077", "47A"
    amap = {k: a for k, a in anchors.load_anchors(v, set()).items() if k != s}
    placed = {k: (a.theta, a.t) for k, a in amap.items()}
    bodies = engine.placed_bodies(v, placed)
    truth = evaluate.truth_pose(v, s)
    assert engine.printed_contradictions(v, s, (truth["theta"], truth["t"]), bodies) == []
    wrong = (truth["theta"], truth["t"] + np.array([129.0, 0.0]))
    assert engine.printed_contradictions(v, s, wrong, bodies), "129 m away reaches nothing"


def test_overlap_ignores_a_strip_as_thin_as_the_anchors_disagree(tmp_path, monkeypatch):
    """47A's correct pose overlapped 171 by a 2 m strip along 120 m; that is not sitting on it."""
    import geopandas as gpd
    from shapely.geometry import Polygon
    from autogeoref import paths
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "V").mkdir(parents=True)
    sheet = Polygon([(0, 0), (120, 0), (120, 25), (0, 25)])
    gpd.GeoDataFrame({"poly_id": [1]}, geometry=[sheet]).to_file(
        tmp_path / "FMB_Vector" / "V" / "A_parcels.geojson", driver="GeoJSON")
    pose = (0.0, np.array([0.0, 0.0]))
    neighbour_2m_over = Polygon([(0, 23), (120, 23), (120, 60), (0, 60)])      # 2 m strip overlap
    neighbour_on_top = Polygon([(10, 5), (110, 5), (110, 20), (10, 20)])       # really on top
    assert engine._overlap("V", "A", pose, {"B": neighbour_2m_over}) <= engine.MAX_OVERLAP_SHARE
    assert engine._overlap("V", "A", pose, {"B": neighbour_on_top}) > engine.MAX_OVERLAP_SHARE


def test_calibration_of_2026_09_19_holds():
    """The 15 Kizhikaranai leave-one-out rows: every green within 2.7 m, the 16 m parcel red."""
    good3 = {"boundary_rms_m": 0.884, "anchor_rms_m": 2.2, "share": 0.667, "observable": False,
             "n_neighbours": 3, "method": "neighbour", "ambiguous": False,
             "support_m": 125.0, "support_next_m": 0.0}
    assert engine.colour_of(good3) == "green"                                   # 43A, 0.58 m
    two = dict(good3, n_neighbours=2, boundary_rms_m=0.165, share=0.482, observable=True)
    assert engine.colour_of(two) == "amber"                                     # 47B, 8.9 m
    one = dict(good3, n_neighbours=1, boundary_rms_m=0.241, share=0.294, observable=True)
    assert engine.colour_of(one) == "red"                                       # 40B, 16.4 m
    imaged_two = dict(two, share=0.811, boundary_rms_m=2.059)
    assert engine.colour_of(imaged_two) == "amber"                              # 48A, 4.1 m


def test_green_needs_a_direct_tie_to_the_teams_own_parcels():
    """47B, seed run 2026-09-19: three neighbours, all placed in the same pass from one wrong guess."""
    row = {"boundary_rms_m": 0.265, "anchor_rms_m": 0.66, "share": 0.5, "observable": False,
           "n_neighbours": 3, "method": "neighbour", "ambiguous": False,
           "support_m": 80.0, "support_next_m": 54.0, "contradictions": 0, "anchor_partners": 0}
    assert engine.colour_of(row) != "green"
    assert engine.colour_of(dict(row, anchor_partners=1)) == "green"


def test_side_checks_read_both_sheets():
    """48A prints 47B to its south and 47B prints 48A to its north; a pose that puts 47B north
    of 48A contradicts both."""
    from autogeoref import anchors, evaluate, paths
    if not paths.vector_dir("35_04_077").exists():
        pytest.skip("village not present")
    v = "35_04_077"
    amap = {k: a for k, a in anchors.load_anchors(v, set()).items() if k == "48A"}
    placed = {k: (a.theta, a.t) for k, a in amap.items()}
    bodies = engine.placed_bodies(v, placed)
    truth = evaluate.truth_pose(v, "47B")
    ok, bad = engine.side_checks(v, "47B", (truth["theta"], truth["t"]), placed, bodies)
    assert ok == 2 and bad == 0
    flipped = (truth["theta"], truth["t"] + np.array([0.0, 160.0]))     # north of 48A instead
    ok2, bad2 = engine.side_checks(v, "47B", flipped, placed, bodies)
    assert bad2 == 2
