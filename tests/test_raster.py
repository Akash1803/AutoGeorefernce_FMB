import pytest
import rasterio

from autogeoref import raster

SPIKE = r"D:\code\FMB_to_GeoJSON\autogeoref\scratch_2026-09-17\spike_gcp\kizhi_sat_z20.tif"


def test_tile_description_is_a_tms_service_for_the_project_layer():
    assert "mt1.google.com/vt/lyrs=s" in raster.TILE_XML
    assert "<TileLevel>20</TileLevel>" in raster.TILE_XML
    assert "EPSG:3857" in raster.TILE_XML


def test_the_spike_raster_is_readable_and_georeferenced():
    with rasterio.open(SPIKE) as src:
        assert src.crs.to_epsg() == 32644
        assert src.res[0] == pytest.approx(0.15, abs=0.01)
        assert src.count == 3


@pytest.mark.network
def test_export_writes_a_small_window_and_pins_it(tmp_path, monkeypatch):
    monkeypatch.setattr(raster.paths, "PROJECT", tmp_path)
    bounds = (392900.0, 1413300.0, 392980.0, 1413380.0)     # 80 x 80 m over Kizhikaranai
    out = raster.export("35_04_077", bounds, zoom=20, resolution=0.15)
    assert out.exists()
    with rasterio.open(out) as src:
        assert src.crs.to_epsg() == 32644
        assert src.width == pytest.approx(533, abs=2) and src.height == pytest.approx(533, abs=2)
        assert src.read(1).any(), "tiles actually downloaded"
    pin = raster.pinned("35_04_077")
    assert pin["path"] == str(out) and len(pin["sha256"]) == 64 and pin["zoom"] == 20


@pytest.mark.network
def test_export_archives_the_previous_raster(tmp_path, monkeypatch):
    monkeypatch.setattr(raster.paths, "PROJECT", tmp_path)
    bounds = (392900.0, 1413300.0, 392940.0, 1413340.0)
    first = raster.export("35_04_077", bounds, run_id="20260101")
    assert first.exists()
    second = raster.export("35_04_077", bounds, run_id="20260102")
    archived = list((tmp_path / "_logs").glob("**/satellite_*.tif"))
    assert archived, "the previous raster is archived, never deleted"
    assert second.exists() and second != first


@pytest.mark.network
def test_open_pinned_refuses_a_raster_that_changed_on_disk(tmp_path, monkeypatch):
    monkeypatch.setattr(raster.paths, "PROJECT", tmp_path)
    out = raster.export("35_04_077", (392900.0, 1413300.0, 392930.0, 1413330.0))
    with open(out, "ab") as fh:
        fh.write(b"tampered")
    with pytest.raises(ValueError):
        raster.open_pinned("35_04_077")
