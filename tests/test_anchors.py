import numpy as np
import pytest

from autogeoref import anchors, paths

VILLAGE = "35_04_077"


def test_read_points_parses_the_qgis_georeferencer_file():
    P = anchors.read_points(paths.points_path(VILLAGE, "47B"))
    assert P.shape[1] == 4 and len(P) >= 3
    assert 380000 < P[0, 0] < 400000 and 1400000 < P[0, 1] < 1420000, "map side is EPSG:32644"
    assert abs(P[0, 2]) < 5000 and abs(P[0, 3]) < 5000, "source side is sheet metres"


def test_write_points_round_trips(tmp_path):
    rows = [(392809.1, 1413357.1, 198.4, 399.9), (392900.0, 1413400.0, 250.0, 420.0)]
    p = tmp_path / "x_parcels.geojson.points"
    anchors.write_points(p, rows, crs_wkt='PROJCRS["WGS 84 / UTM zone 44N"]')
    text = p.read_text(encoding="utf-8").splitlines()
    assert text[0].startswith("#CRS: ")
    assert text[1] == "mapX,mapY,sourceX,sourceY,enable,dX,dY,residual"
    back = anchors.read_points(p)
    assert np.allclose(back, np.array(rows), atol=1e-6)


def test_pose_from_points_is_rigid_and_matches_the_hand_file_roughly():
    got = anchors.pose_from_points(VILLAGE, "46A")
    assert got is not None
    theta, t, rms, scale, aniso, n = got
    assert n >= 3
    assert rms < 1.0, "46A was placed cleanly; its own GCPs fit a rigid pose"
    assert 0.8 < scale < 1.2


def test_load_anchors_takes_one_file_per_survey_and_flags_the_stretch():
    got = anchors.load_anchors(VILLAGE, tool_written=set())
    assert len(got) == 15, "all 15 Kizhikaranai surveys are hand placed"
    assert got["47A"].file.name == "47A_parcels_modified2.gpkg", "newest variant wins"
    assert got["48B"].scale < 0.95, "48B is georeferenced ~10 % small"
    assert all(a.source in ("points", "geometry") for a in got.values())


def test_fingerprint_is_stable_and_geometry_sensitive(tmp_path):
    import geopandas as gpd
    from shapely.geometry import Polygon
    a = tmp_path / "a.gpkg"
    gpd.GeoDataFrame({"x": [1]}, geometry=[Polygon([(0, 0), (1, 0), (1, 1)])], crs="EPSG:32644").to_file(a, driver="GPKG")
    f1 = anchors.fingerprint(a)
    assert f1 == anchors.fingerprint(a)
    b = tmp_path / "b.gpkg"
    gpd.GeoDataFrame({"x": [1]}, geometry=[Polygon([(0, 0), (1, 0), (1, 1.5)])], crs="EPSG:32644").to_file(b, driver="GPKG")
    assert anchors.fingerprint(b) != f1


def test_anchor_conflicts_only_reports_pairs_that_exceed_their_own_tolerance():
    got = anchors.load_anchors(VILLAGE, tool_written=set())
    conf = anchors.anchor_conflicts(VILLAGE, got)
    assert conf, "Kizhikaranai has anchors that disagree; the check must surface them"
    for c in conf:
        assert set(c) >= {"a", "b", "p90_gap_m", "allowed_m", "overlap_sqm"}
        assert c["p90_gap_m"] > c["allowed_m"], "a pair inside its tolerance is not a conflict"
    pairs = {tuple(sorted((c["a"], c["b"]))) for c in conf}
    assert ("170", "42A") in pairs, "the widest disagreement on this village"
    assert ("43A", "46B") not in pairs, "these two share an exact 27.4 m boundary"


def test_48A_and_171_agree_once_the_pose_comes_from_the_teams_own_gcps():
    """The 2026-09-17 prototype put these 10.7 m apart; that was its cruder pose, not the data."""
    got = anchors.load_anchors(VILLAGE, tool_written=set())
    a = anchors.placed_geometry(VILLAGE, got["48A"])
    b = anchors.placed_geometry(VILLAGE, got["171"])
    assert a.distance(b) == 0.0
    assert a.intersection(b).area < 100.0
