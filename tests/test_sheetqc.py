import json

from shapely.geometry import Polygon

from autogeoref import paths, sheetqc


def _sq(x0, y0, x1, y1):
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def test_check_plots_flags_sliver_spike_unnumbered_and_overlap():
    plots = [("1", _sq(0, 0, 20, 20)), ("None", _sq(20, 0, 20.1, 20)),
             ("2", Polygon([(30, 0), (60, 0), (30, 2)])), ("3", _sq(10, 10, 25, 25))]
    issues = sheetqc.check_plots(plots)
    text = " | ".join(issues)
    assert "sliver" in text and "unnumbered plot" in text
    assert "spike" in text
    assert "plots 1 and 3 overlap" in text


def test_check_plots_is_quiet_on_a_clean_sheet():
    assert sheetqc.check_plots([("1", _sq(0, 0, 20, 20)), ("2", _sq(20, 0, 40, 20))]) == []


def test_check_lines_flags_a_thin_line_crossing_the_boundary():
    def line(layer, a, b):
        return {"properties": {"layer": layer}, "geometry": {"type": "LineString", "coordinates": [a, b]}}
    lines = [line("survey_boundary_main", (0, 0), (0, 100)), line("subdivision_line", (5, 50), (-5, 60)),
             line("subdivision_line", (0, 20), (30, 20))]
    issues = sheetqc.check_lines(lines)
    assert len(issues) == 1 and "crosses the boundary" in issues[0]


def test_qc_village_writes_a_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    d = tmp_path / "FMB_Vector" / "v"
    d.mkdir(parents=True)
    feats = [{"type": "Feature", "properties": {"plot_no": "1"}, "geometry": _sq(0, 0, 20, 20).__geo_interface__},
             {"type": "Feature", "properties": {"plot_no": None}, "geometry": _sq(20, 0, 20.2, 20).__geo_interface__}]
    (d / "7_parcels.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")
    res = sheetqc.qc_village("v")
    assert res["sheets"][0]["survey"] == "7" and res["sheets"][0]["n_issues"] >= 2
    assert (d / "sheet_qc.csv").read_text(encoding="utf-8").count("\n") == 2
