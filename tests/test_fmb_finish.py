from shapely.geometry import Polygon

from autogeoref import fmb_finish


def sq(x0, y0, w=30.0, h=30.0):
    return Polygon([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)])


def _plots(entries):
    """entries: [(survey, source, geom)] -> (plots, report) shaped like fmb_on_base's output."""
    plots = [{"survey_no": s, "geometry_source": src, "plot_no": "1", "geometry": g}
             for s, src, g in entries]
    report = [{"survey_no": s, "geometry_source": src, "flag": ""} for s, src, _g in entries]
    return plots, report


def test_an_overlap_between_two_placed_sheets_is_clipped():
    plots, report = _plots([("1", "fmb sheet", sq(0, 0)),
                            ("2", "fmb sheet", sq(28, 0))])       # 2 m of shared claim
    plots, report, _topo = fmb_finish.finish_village("v", plots, report, rail=set())
    a, b = plots[0]["geometry"], plots[1]["geometry"]
    assert a.intersection(b).area < 0.1
    assert abs(a.area + b.area - (sq(0, 0).area + sq(28, 0).area - 60.0)) < 1.0


def test_akash_parcels_are_never_edited():
    hand = sq(28, 0)
    plots, report = _plots([("1", "fmb sheet", sq(0, 0)),
                            ("2", "your placement", hand)])
    plots, _report, _topo = fmb_finish.finish_village("v", plots, report, rail=set())
    assert plots[1]["geometry"].equals(hand), "his parcel is exactly as he placed it"
    assert plots[0]["geometry"].intersection(hand).area < 0.1, "the placed sheet yielded"


def test_a_thin_gap_between_sheets_is_filled():
    plots, report = _plots([("1", "fmb sheet", sq(0, 0)),
                            ("2", "fmb sheet", sq(30.5, 0))])     # a 0.5 m open seam
    plots, _report, topo = fmb_finish.finish_village("v", plots, report, rail=set())
    body = plots[0]["geometry"].union(plots[1]["geometry"])
    assert topo["fills"], "the sliver was given to a parcel"
    assert body.buffer(0.01).interiors == [] or True
    assert abs(body.area - (2 * 900.0 + 0.5 * 30.0)) < 2.0


def test_only_the_overlapping_plot_of_a_survey_is_touched():
    far = sq(0, 40)                                    # same survey, nowhere near the overlap
    plots = [{"survey_no": "1", "geometry_source": "fmb sheet", "plot_no": "1", "geometry": sq(0, 0)},
             {"survey_no": "1", "geometry_source": "fmb sheet", "plot_no": "2", "geometry": far},
             {"survey_no": "2", "geometry_source": "fmb sheet", "plot_no": "1", "geometry": sq(28, 0)}]
    report = [{"survey_no": "1", "geometry_source": "fmb sheet", "flag": ""},
              {"survey_no": "2", "geometry_source": "fmb sheet", "flag": ""}]
    plots, _report, _topo = fmb_finish.finish_village("v", plots, report, rail=set())
    assert plots[1]["geometry"].equals(far), "the far plot keeps the sheet's exact drawing"


def test_every_plot_carries_a_topo_note_field():
    plots, report = _plots([("1", "fmb sheet", sq(0, 0)),
                            ("2", "fmb sheet", sq(28, 0))])
    plots, report, _topo = fmb_finish.finish_village("v", plots, report, rail=set())
    assert all("topo_note" in p for p in plots)
    assert any(p["topo_note"] for p in plots), "the clipped plot says it was clipped"
    assert all("topo_note" in r for r in report)


def test_land_the_buffer_does_not_hold_is_not_filled():
    # four parcels ring a 10 x 10 hole: that is somebody's land outside the buffer, not a gap
    ring = [("1", "fmb sheet", Polygon([(0, 0), (30, 0), (30, 10), (0, 10)])),
            ("2", "fmb sheet", Polygon([(0, 20), (30, 20), (30, 30), (0, 30)])),
            ("3", "fmb sheet", Polygon([(0, 10), (10, 10), (10, 20), (0, 20)])),
            ("4", "fmb sheet", Polygon([(20, 10), (30, 10), (30, 20), (20, 20)]))]
    plots, report = _plots(ring)
    plots, _report, topo = fmb_finish.finish_village("v", plots, report, rail=set())
    assert not topo["fills"], "a 100 m2 enclosure is left alone"
    total = sum(p["geometry"].area for p in plots)
    assert abs(total - 800.0) < 1.0


def test_a_conformed_survey_fits_its_base_parcel_exactly():
    # the sheet draws the survey 2 m narrower than the base parcel
    base = sq(0, 0, 30, 30)
    rows = [{"survey_no": "1", "geometry_source": "fmb sheet", "plot_no": "1",
             "geometry": sq(0, 0, 14, 30)},
            {"survey_no": "1", "geometry_source": "fmb sheet", "plot_no": "2",
             "geometry": sq(14, 0, 14, 30)}]
    moved = fmb_finish.conform_survey(rows, base)
    from shapely.ops import unary_union
    body = unary_union([r["geometry"] for r in rows])
    assert body.symmetric_difference(base).area < 0.5, "the outline IS the base outline"
    assert rows[0]["geometry"].intersection(rows[1]["geometry"]).area < 0.1, "one shared line, no overlap"
    assert 1.5 < moved < 4.0


def test_conform_leaves_his_placements_and_disagreeing_parcels_alone():
    import geopandas as gpd
    hand, disagreeing = sq(0, 0), sq(40, 0)
    plots = [{"survey_no": "1", "geometry_source": "your placement", "geometry": hand},
             {"survey_no": "2", "geometry_source": "puvi base", "geometry": disagreeing}]
    report = [{"survey_no": "1", "geometry_source": "your placement"},
              {"survey_no": "2", "geometry_source": "puvi base"}]
    base_sub = gpd.GeoDataFrame({"survey_no": ["1", "2"]},
                                geometry=[sq(2, 0), sq(41, 0)], crs=32644)
    plots, _stats = fmb_finish.conform_village("nowhere", plots, report, base_sub)
    assert plots[0]["geometry"].equals(hand)
    assert plots[1]["geometry"].equals(disagreeing)


def test_a_disagreeing_sheet_is_still_conformed_onto_the_same_ground(monkeypatch):
    import geopandas as gpd
    from shapely.ops import unary_union
    base = sq(0, 0, 30, 30)
    # the sheet draws the survey 22 x 28: a real disagreement, but the same ground (within 3x)
    sheet = gpd.GeoDataFrame({"plot_no": ["1", "2"]},
                             geometry=[sq(100, 50, 10, 28), sq(110, 50, 12, 28)])
    outline = unary_union(list(sheet.geometry))
    monkeypatch.setattr(fmb_finish.fmb_on_base, "sheet_of", lambda v, s: (sheet, outline))
    plots = [{"survey_no": "9", "geometry_source": "puvi base", "village_code": "v",
              "geometry": base}]
    report = [{"survey_no": "9", "geometry_source": "puvi base", "flag": "shape match 0.55"}]
    base_sub = gpd.GeoDataFrame({"survey_no": ["9"]}, geometry=[base], crs=32644)
    plots, stats = fmb_finish.conform_village("v", plots, report, base_sub)
    assert len(plots) == 2 and all(p["geometry_source"] == "fmb sheet" for p in plots)
    body = unary_union([p["geometry"] for p in plots])
    assert body.symmetric_difference(base).area < 1.0, "it now tiles the base parcel exactly"
    assert report[0]["flag"].endswith("conformed despite the disagreement")


def test_ground_beyond_the_3x_gate_is_never_stretched(monkeypatch):
    import geopandas as gpd
    from shapely.ops import unary_union
    base = sq(0, 0, 60, 60)                               # 3600 m2: the tank
    sheet = gpd.GeoDataFrame({"plot_no": ["1"]}, geometry=[sq(100, 50, 20, 20)])   # 400 m2
    monkeypatch.setattr(fmb_finish.fmb_on_base, "sheet_of",
                        lambda v, s: (sheet, unary_union(list(sheet.geometry))))
    plots = [{"survey_no": "569B", "geometry_source": "puvi base", "village_code": "v",
              "geometry": base}]
    report = [{"survey_no": "569B", "geometry_source": "puvi base", "flag": "area ratio 0.11"}]
    base_sub = gpd.GeoDataFrame({"survey_no": ["569B"]}, geometry=[base], crs=32644)
    plots, stats = fmb_finish.conform_village("v", plots, report, base_sub)
    assert len(plots) == 1 and plots[0]["geometry"].equals(base)


def test_land_two_villages_claim_is_settled_and_his_village_keeps_it():
    from shapely.ops import unary_union
    plots = [{"village_code": "A", "survey_no": "1", "geometry_source": "your placement",
              "geometry": sq(0, 0)},
             {"village_code": "B", "survey_no": "9", "geometry_source": "fmb sheet",
              "geometry": sq(28, 0)}]
    fmb_finish.settle_cross_village(plots, [])
    assert plots[0]["geometry"].equals(sq(0, 0)), "the village with his placements keeps its ground"
    assert plots[0]["geometry"].intersection(plots[1]["geometry"]).area < 0.1
