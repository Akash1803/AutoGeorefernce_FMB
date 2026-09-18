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
