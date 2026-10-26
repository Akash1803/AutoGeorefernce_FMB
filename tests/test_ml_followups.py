"""Akash's follow-ups of 2026-09-21 on PR 0 and PR 1."""
import csv
import json

from autogeoref import evalrows, paths, report, transcribe
from autogeoref.transcribe import Entry


def _row(village, survey, run_id, label, run_at="2026-09-21T10:00:00", **kw):
    d = {"run_id": run_id, "run_at": run_at, "village": village, "survey": survey, "colour": "green",
         "err_fmb_m": 1.0, "err_hand_m": 1.0, "within_3m": label, "stretched": "0", "notes": "", "sheet_qc_reason": ""}
    d.update(kw)
    return d


def test_distinct_summary_counts_parcels_not_rows_and_classes_by_latest_row():
    rows = [_row("A", "1", "r1", "1"), _row("A", "1", "r2", "0", run_at="2026-09-21T11:00:00"),
            _row("A", "2", "r1", "1"), _row("A", "3", "r1", ""), _row("B", "1", "r1", "0")]
    d = evalrows.distinct_summary(rows)
    assert d["rows"] == 5 and d["labelled_rows"] == 4
    assert d["distinct_parcels"] == 4 and d["distinct_labelled"] == 3
    assert d["parcels_within"] == 1 and d["parcels_beyond"] == 2, "parcel A/1 counts once, by its latest row"
    assert d["villages"] == ["A", "B"]


def test_pr3_gate_is_on_distinct_parcels_across_two_villages():
    rows = [_row("A", str(i), "r1", "1") for i in range(300)] + [_row("A", str(i), "r2", "1") for i in range(300)]
    g = report.pr3_gate(rows)
    assert g["distinct_labelled"] == 300 and not g["ok"], "600 rows, 300 parcels, one village"
    rows.append(_row("B", "1", "r1", "0"))
    assert report.pr3_gate(rows)["ok"]


def test_append_rows_dedupes_on_run_id_and_migrates_the_header(tmp_path):
    p = tmp_path / "rows.csv"
    old_fields = [f for f in evalrows.FIELDS if f != "sheet_qc_reason"]
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=old_fields)
        w.writeheader()
        w.writerow({k: "" for k in old_fields} | {"run_id": "r1", "village": "A", "survey": "1"})
    dup = evalrows.qc_failure_row("A", "1", "x", "r1")
    new = evalrows.qc_failure_row("A", "2", "drawn 9.9 m short", "r1")
    evalrows.append_rows([dup, new], p)
    rows = evalrows.load_rows(p)
    assert [r["survey"] for r in rows] == ["1", "2"], "the duplicate of r1/A/1 was skipped"
    assert list(rows[0].keys()) == evalrows.FIELDS, "header migrated to the current schema"
    assert rows[1]["sheet_qc_reason"] == "drawn 9.9 m short" and rows[1]["within_3m"] == ""


def test_qc_failure_rows_carry_no_label_and_stay_out_of_the_statistics():
    rows = [_row("A", "1", "r1", "1"), _row("A", "2", "r1", "0", err_fmb_m=40.0),
            evalrows.qc_failure_row("A", "613", "east side drawn 60.9 m, printed 70.8 m", "qc").as_dict()]
    assert not evalrows.is_labelled(rows[2])
    m = report.metrics(rows)
    assert m["rows"] == 2 and m["labelled"] == 2 and m["distinct_labelled"] == 2
    assert all(w["survey"] != "613" for w in report.worst10(rows))
    text = report.write_report(rows, __import__("pathlib").Path(__import__("tempfile").mkdtemp()) / "r.md", "t").read_text(encoding="utf-8")
    assert "Sheet-QC failures" in text and "A 613: east side" in text and "Gate (300 parcels, 2 villages): not met" in text


def test_transcribed_fields_carry_a_source_and_resolutions_are_hand():
    a = {"91": [Entry("7", "SE"), Entry("44", "S")]}
    b = {"91": [Entry("44", "SW")]}
    table, dis, _ = transcribe.merge(a, b, ["91"], "v")
    assert all(r["source"] == "reader" for r in table["91"])
    for d in dis:
        d.resolution = "yes" if d.field == "number" else "S"
    transcribe.apply_resolutions(table, dis)
    seven = next(r for r in table["91"] if r["number"] == "7")
    forty = next(r for r in table["91"] if r["number"] == "44")
    assert seven["source"] == "hand" and seven["resolved_by"] == "hand" and seven["readers"] == 2
    assert forty["source"] == "hand" and forty["side"] == "S"
    assert all(d.resolved_by == "hand" for d in dis)


def test_review_csv_has_resolved_by_column(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "v").mkdir(parents=True)
    a = {"91": [Entry("7", "SE")]}
    _, dis, _ = transcribe.merge(a, {"91": []}, ["91"], "v")
    dis[0].resolution = "no"
    p = transcribe.write_review_csv("v", dis)
    rows = list(csv.DictReader(p.open(encoding="utf-8")))
    assert rows[0]["resolved_by"] == "hand" and "resolved_by" in transcribe.REVIEW_COLUMNS


def test_hand_corrected_sheets_are_excluded_from_seeds(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    d = tmp_path / "FMB_Vector" / "v"
    d.mkdir(parents=True)
    (d / "sheet_corrections.csv").write_text(
        "date,survey,edge,drawn_m,printed_m,also_printed_on,change,area_before_m2,area_after_m2,backup,verified_by\n"
        "x,613,east side,60.9,70.8,569C,moved,1,2,b,eye\nx,47A,plots,,,,re-converted,1,1,b,converter\n", encoding="utf-8")
    assert transcribe.corrected_sheets("v") == {"613"}
    table = {"613": [{"number": "614", "side": "S", "readers": 2, "confidence": "high"},
                     {"number": "607", "side": "N", "readers": 2, "confidence": "high"}],
             "614": [], "607": []}
    plan = transcribe.seed_plan("v", table, ["613", "614", "607"], placed=["613", "614"])
    assert "613" not in plan["from_scratch"]["seeds"] and "613" not in plan["additional"]["seeds"]
    assert plan["excluded_from_seeds_and_labels"] == ["613"]
    assert "613" not in plan["from_scratch"]["labelled_parcels"] and "613" not in plan["additional"]["labelled_parcels"]
    graph = {"613": {"614", "607"}, "614": {"613"}, "607": {"613"}}
    seeds = transcribe.dominating_set(graph, exclude=frozenset({"613"}))
    assert "613" not in seeds and set(seeds) == {"614", "607"}
