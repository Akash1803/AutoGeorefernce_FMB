import numpy as np
import pytest

from autogeoref import engine


def test_colour_rule_green_needs_an_anchor_and_observability():
    assert engine.colour_of({"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.65,
                             "observable": True, "n_neighbours": 2, "method": "anchors",
                             "ambiguous": False}) == "green"


def test_colour_rule_amber_when_only_one_condition_is_short():
    assert engine.colour_of({"boundary_rms_m": 0.20, "anchor_rms_m": 0.3, "share": 0.20,
                             "observable": False, "n_neighbours": 1, "method": "anchors",
                             "ambiguous": False}) == "amber"


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


def test_run_on_a_village_where_every_parcel_is_an_anchor_changes_nothing():
    """Kizhikaranai is fully hand placed: the engine must place nothing and touch no manual file."""
    from autogeoref import paths
    if not paths.vector_dir("35_04_077").exists():
        pytest.skip("village not present")
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
