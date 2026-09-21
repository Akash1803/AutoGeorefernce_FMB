import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from autogeoref import files, paths, qgis_bridge

SQ = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])


def _write(path, value):
    gpd.GeoDataFrame({"v": [value]}, geometry=[SQ], crs="EPSG:32644").to_file(path, driver="GPKG")


def test_safe_write_replaces_a_file_that_nothing_holds(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)
    assert files.safe_write(target, lambda p: _write(p, 2))
    assert int(gpd.read_file(target)["v"].iloc[0]) == 2


def test_safe_write_leaves_the_original_when_the_writer_raises(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)

    def boom(p):
        raise RuntimeError("writer failed")

    assert files.safe_write(target, boom) is False
    assert int(gpd.read_file(target)["v"].iloc[0]) == 1
    assert not list(tmp_path.glob("*autogeoref-tmp*"))


def test_safe_write_clears_stale_wal_sidecars(tmp_path):
    target = tmp_path / "x.gpkg"
    _write(target, 1)
    (tmp_path / "x.gpkg-wal").write_bytes(b"stale")
    (tmp_path / "x.gpkg-shm").write_bytes(b"stale")
    assert files.safe_write(target, lambda p: _write(p, 2))
    assert not (tmp_path / "x.gpkg-wal").exists()


def test_supersede_moves_and_never_deletes(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    target = tmp_path / "old.gpkg"
    _write(target, 1)
    moved = files.supersede(target, "20260918")
    assert moved is not None and moved.exists() and not target.exists()
    assert "georef_superseded_20260918" in str(moved)


def test_held_layers_is_empty_when_qgis_is_unreachable(monkeypatch, tmp_path):
    monkeypatch.setattr(qgis_bridge, "available", lambda: False)
    assert files.held_layers(tmp_path / "anything.gpkg") == []


@pytest.mark.skipif(not qgis_bridge.available(), reason="QGIS is not running")
def test_held_layers_finds_a_file_the_open_project_holds():
    p = paths.vector_dir("35_04_077") / "171_parcels_modified1.gpkg"
    got = files.held_layers(p)
    assert isinstance(got, list)
    for layer in got:
        assert {"id", "name", "editable"} <= set(layer)


@pytest.mark.skipif(not qgis_bridge.available(), reason="QGIS is not running")
def test_the_bridge_runs_code_inside_qgis():
    out = qgis_bridge.json_result("import json; print(json.dumps({'ok': 1 + 1}))")
    assert out == {"ok": 2}



def test_safe_write_refuses_to_take_a_held_layer_out_of_the_project_by_default(tmp_path, monkeypatch):
    """2026-09-21: replacing ten of the team's files removed their layers from the open project."""
    target = tmp_path / "x.gpkg"
    _write(target, 1)
    monkeypatch.setattr(files, "held_layers", lambda p: [{"id": "L1", "name": "x", "editable": False}])
    released = []
    monkeypatch.setattr(files, "release", lambda p: released.append(p) or [])
    assert files.safe_write(target, lambda p: _write(p, 2)) is False
    assert released == [], "no layer was removed"
    assert int(gpd.read_file(target)["v"].iloc[0]) == 1, "the file is untouched"
    assert files.safe_write(target, lambda p: _write(p, 2), allow_release=True) is True
