import csv
import math

import geopandas as gpd
import pytest
from shapely import affinity
from shapely.geometry import Polygon

from autogeoref import config, paths, references

V, S = "99_99_901", "7A"
THETA, T = 30.0, (500000.0, 1400000.0)


def _sheet():
    return [({"poly_id": 1, "plot_no": "1"}, Polygon([(0, 0), (60, 0), (60, 40), (0, 40)])),
            ({"poly_id": 2, "plot_no": "2"}, Polygon([(60, 0), (100, 0), (100, 40), (60, 40)]))]


def _place(geom, sx=1.0):
    g = affinity.scale(geom, xfact=sx, yfact=1.0, origin=(0, 0))
    g = affinity.rotate(g, THETA, origin=(0, 0))
    return affinity.translate(g, T[0], T[1])


@pytest.fixture
def village(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    d = tmp_path / "FMB_Vector" / V
    d.mkdir(parents=True)
    sheet = _sheet()
    gpd.GeoDataFrame([p for p, _g in sheet], geometry=[g for _p, g in sheet]).to_file(
        d / ("%s_parcels.geojson" % S), driver="GeoJSON")
    return d


def _hand(d, sx):
    sheet = _sheet()
    gpd.GeoDataFrame([p for p, _g in sheet], geometry=[_place(g, sx) for _p, g in sheet],
                     crs="EPSG:32644").to_file(d / ("%s_parcels_modified.gpkg" % S), driver="GPKG")


def test_an_unstretched_hand_placement_is_a_clean_reference(village):
    _hand(village, 1.0)
    cfg = config.Config()
    ref = references.reference_for(V, S, cfg)
    assert ref is not None and ref.source == "geometry"
    assert ref.area_scale == pytest.approx(1.0, abs=1e-3)
    assert ref.disagree_mean_m < 0.05 and ref.disagree_max_m < 0.1
    assert ref.disagree_corner_mean_m < 0.05 and ref.corners_matched >= 4
    assert abs(ref.stretch_long_pct) < 0.5 and abs(ref.stretch_short_pct) < 0.5
    assert ref.band == "<=5%" and not ref.stretched and not ref.disputed
    assert math.isclose(ref.geom_fmb.area, 4000.0, rel_tol=1e-6)


def test_a_stretched_hand_placement_is_excluded_by_every_rule_and_kept_as_a_feature(village):
    _hand(village, 1.25)                                   # 25 % longer along the sheet's x axis
    for rule in ("pct", "axis", "displacement"):
        cfg = config.Config()
        cfg.stretch.rule = rule
        ref = references.reference_for(V, S, cfg)
        assert ref.stretched, rule
        assert ref.stretched_reason
    assert ref.stretch_long_pct == pytest.approx(25.0, abs=1.0)
    assert abs(ref.stretch_short_pct) < 1.0
    assert ref.area_scale == pytest.approx(math.sqrt(1.25), abs=0.01)
    assert ref.disagree_mean_m < ref.disagree_corner_mean_m, "a slide along the long edges hides from outline sampling"
    assert ref.disagree_corner_mean_m > 3.0, "the corners of a 100 m parcel stretched by 25 % move by metres"


def test_the_middle_tier_is_labelled_under_the_displacement_rule_when_it_moves_little(village):
    """A 6 % stretch of a 100 m parcel: the corners move under 3 m on average, so it keeps its label."""
    _hand(village, 1.06)
    cfg = config.Config()
    cfg.stretch.rule = "displacement"
    ref = references.reference_for(V, S, cfg)
    assert ref.band in ("<=5%", "5-10%")
    assert not ref.stretched
    cfg.stretch.rule = "axis"
    assert not references.reference_for(V, S, cfg).stretched, "6 % is under the 10 % limit"


def test_disputed_surveys_come_from_the_village_list_and_the_audit(village):
    _hand(village, 1.0)
    (village / "disputed.csv").write_text("survey,reason\n7A,team says wrong sheet\n", encoding="utf-8")
    ref = references.reference_for(V, S, config.Config())
    assert ref.disputed and "wrong sheet" in ref.disputed_reason


def test_write_exclusions_produces_the_three_corridor_files_and_keeps_other_villages(village, tmp_path):
    _hand(village, 1.25)
    cfg = config.Config()
    refs = references.references(V, cfg)
    assert refs[S].stretched
    out = tmp_path / "eval"
    out.mkdir()
    other = {c: "" for c in references.STRETCHED_COLUMNS}
    other.update(village="other", survey="1", rule="pct", reason="x")
    with (out / "excluded_stretched.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=references.STRETCHED_COLUMNS)
        w.writeheader()
        w.writerow(other)
    p1, p2, p3 = references.write_exclusions(V, refs, cfg, out)
    rows = list(csv.DictReader(p1.open(encoding="utf-8")))
    assert {r["village"] for r in rows} == {"other", V}, "another village's rows survive"
    assert any(r["survey"] == S and r["rule"] == cfg.stretch.rule for r in rows)
    assert list(csv.DictReader(p2.open(encoding="utf-8"))) == []
    ref_rows = list(csv.DictReader(p3.open(encoding="utf-8")))
    assert ref_rows[0]["survey"] == S and ref_rows[0]["stretched"] == "True"
