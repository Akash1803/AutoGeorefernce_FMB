import csv
import json

import pytest

from autogeoref import config, paths, transcribe
from autogeoref.transcribe import Entry


def _reading(**sheets):
    return {s: [Entry(*e) for e in entries] for s, entries in sheets.items()}


def test_agreed_numbers_are_trusted_and_slash_forms_match():
    a = _reading(**{"584": [("569/C", "NE"), ("614", "N"), ("101A", "NE")]})
    b = _reading(**{"584": [("569C", "NE"), ("614", "N"), ("101A", "NE")]})
    table, dis, stats = transcribe.merge(a, b, ["584"], "v")
    assert {(r["number"], r["side"], r["readers"]) for r in table["584"]} == {("569C", "NE", 2), ("614", "N", 2), ("101A", "NE", 2)}
    assert dis == []
    assert stats[0].n_both == 3 and stats[0].n_union == 3 and stats[0].jaccard == 1.0


def test_side_conflict_keeps_the_number_but_withholds_the_side():
    a = _reading(**{"85A": [("55", "S")]})
    b = _reading(**{"85A": [("55", "SW")]})
    table, dis, _ = transcribe.merge(a, b, ["85A"], "v")
    assert table["85A"] == [{"number": "55", "side": "", "readers": 2, "confidence": "high"}]
    assert len(dis) == 1 and dis[0].field == "side" and dis[0].side_a == "S" and dis[0].side_b == "SW"


def test_number_read_by_one_reader_is_untrusted_and_listed():
    a = _reading(**{"91": [("584", "NW"), ("7", "SE")]})
    b = _reading(**{"91": [("584", "NW")]})
    table, dis, stats = transcribe.merge(a, b, ["91"], "v", "json:A", "json:B")
    seven = next(r for r in table["91"] if r["number"] == "7")
    assert seven["readers"] == 1 and seven["confidence"] == "low"
    assert [d.field for d in dis] == ["number"] and dis[0].note == "read by json:A only"
    assert stats[0].n_a == 2 and stats[0].n_b == 1 and stats[0].n_both == 1


def test_village_labels_and_own_number_are_dropped():
    a = _reading(**{"171": [("V.No. 74 Thirukachur", "E"), ("171", "N"), ("170", "S")]})
    b = _reading(**{"171": [("170", "S")]})
    table, dis, _ = transcribe.merge(a, b, ["171"], "v")
    assert [r["number"] for r in table["171"]] == ["170"] and dis == []


def test_resolutions_fold_back_into_the_table():
    a = _reading(**{"91": [("7", "SE"), ("44", "S")], "85A": [("55", "S")]})
    b = _reading(**{"91": [("44", "SW")], "85A": [("55", "SW")]})
    table, dis, _ = transcribe.merge(a, b, ["91", "85A"], "v")
    for d in dis:
        if d.field == "number":
            d.resolution = "yes"
        elif d.survey == "91":
            d.resolution = "S"
        else:
            d.resolution = "no such side"
    assert transcribe.apply_resolutions(table, dis) == 2
    seven = next(r for r in table["91"] if r["number"] == "7")
    assert seven["readers"] == 2 and seven["verified_by"] == "analyst"
    assert next(r for r in table["91"] if r["number"] == "44")["side"] == "S"
    assert next(r for r in table["85A"] if r["number"] == "55")["side"] == ""


def test_review_csv_round_trip_keeps_resolutions(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "v").mkdir(parents=True)
    a = _reading(**{"91": [("7", "SE")]})
    b = _reading(**{"91": []})
    _, dis, _ = transcribe.merge(a, b, ["91"], "v")
    p = transcribe.write_review_csv("v", dis)
    rows = list(csv.DictReader(p.open(encoding="utf-8")))
    assert rows[0]["field"] == "number" and rows[0]["resolution"] == ""
    rows[0]["resolution"] = "no"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=transcribe.REVIEW_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    table, dis2, _ = transcribe.merge(a, b, ["91"], "v")
    assert transcribe.carry_resolutions(dis2, p) == 1 and dis2[0].resolution == "no"
    transcribe.apply_resolutions(table, dis2)
    assert table["91"] == []


def test_dominating_set_covers_every_node_within_one_hop():
    graph = {"1": {"2"}, "2": {"1", "3"}, "3": {"2", "4"}, "4": {"3", "5"}, "5": {"4"}, "9": set()}
    seeds = transcribe.dominating_set(graph)
    for n in graph:
        assert n in seeds or graph[n] & set(seeds), n
    assert "9" in seeds, "an isolated sheet is its own seed"
    with_fixed = transcribe.dominating_set(graph, fixed=frozenset({"2"}))
    assert "2" not in with_fixed and "1" not in with_fixed
    for n in graph:
        assert n in with_fixed or n == "2" or graph[n] & (set(with_fixed) | {"2"}), n
    assert transcribe.dominating_set(graph, prefer=frozenset({"4"}))[0] in ("2", "4")


def test_neighbour_graph_uses_trusted_entries_only_and_lists_outsiders():
    table = {"584": [{"number": "92A", "side": "E", "readers": 2, "confidence": "high"},
                     {"number": "563", "side": "E", "readers": 2, "confidence": "high"},
                     {"number": "7", "side": "S", "readers": 1, "confidence": "low"}],
             "92A": []}
    graph, outside = transcribe.neighbour_graph(table, ["584", "92A", "7"])
    assert graph["584"] == {"92A"} and graph["92A"] == {"584"} and graph["7"] == set()
    assert outside == {"584": {"563"}}


def test_seed_plan_counts_labelled_parcels_as_placed_minus_seeds():
    table = {"1": [{"number": "2", "side": "E", "readers": 2, "confidence": "high"}],
             "2": [{"number": "3", "side": "E", "readers": 2, "confidence": "high"}],
             "3": []}
    plan = transcribe.seed_plan("v", table, ["1", "2", "3"], placed=["1", "2"])
    assert plan["from_scratch"]["seeds"] == ["2"] and plan["from_scratch"]["labelled_parcels"] == ["1"]
    assert plan["additional"]["seeds"] == [] and plan["edges"] == 2


def test_parse_reading_rejects_a_bad_side():
    with pytest.raises(ValueError, match="side"):
        transcribe.parse_reading({"1": [{"number": "2", "side": "NORTH"}]})
    r = transcribe.parse_reading({"1": [{"number": "2"}, {"number": ""}]})
    assert r["1"] == [Entry("2", "", "high")]


def test_reader_spec_must_be_json_in_pr1():
    assert transcribe.reader_for("json:A").name == "A"
    with pytest.raises(ValueError):
        transcribe.reader_for("paddleocr")


def test_render_sheet_crops_to_the_drawing(tmp_path):
    pymupdf = pytest.importorskip("pymupdf")
    doc = pymupdf.open()
    page = doc.new_page(width=600, height=800)
    page.draw_rect(pymupdf.Rect(200, 300, 400, 500), width=2)
    page.insert_text((250, 290), "101A", fontsize=12)
    pdf = tmp_path / "s.pdf"
    doc.save(str(pdf))
    out = transcribe.render_sheet(pdf, tmp_path / "s.png", scale=1.0, pad_px=20)
    from PIL import Image
    im = Image.open(out)
    assert 200 < im.width < 600 and 200 < im.height < 800, im.size
    assert transcribe.render_sheet(pdf, out, scale=1.0) == out


def test_config_transcription_defaults_and_validation(tmp_path):
    cfg = config.load()
    assert list(cfg.transcription.readers) == ["json:A", "json:B"]
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"transcription": {"readers": ["json:A"]}}), encoding="utf-8")
    with pytest.raises(ValueError, match="two"):
        config.load(str(p))
    p.write_text(json.dumps({"transcription": {"readers": ["json:A", "ocr:B"]}}), encoding="utf-8")
    with pytest.raises(ValueError, match="json"):
        config.load(str(p))
