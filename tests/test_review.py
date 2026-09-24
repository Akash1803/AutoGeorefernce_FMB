import csv
import zipfile

import pytest

from autogeoref import paths, review

ROWS = [{"survey": "47B", "status": "placed", "method": "anchors", "colour": "green",
         "confidence": 88, "share": 0.62, "margin": 0.20, "boundary_rms_m": 0.21,
         "puvi_reference_m": 2.3, "notes": "", "fp_placed": "aa" * 32, "fp_final": "bb" * 32},
        {"survey": "46B", "status": "review", "method": "puvi-only", "colour": "red",
         "confidence": 30, "share": 0.09, "margin": 0.01, "boundary_rms_m": "",
         "puvi_reference_m": 8.1, "notes": "ambiguous", "fp_placed": "", "fp_final": ""}]


def test_status_csv_round_trips(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "35_04_077").mkdir(parents=True)
    p = review.write_status("35_04_077", ROWS)
    assert p.exists()
    back = review.read_status("35_04_077")
    assert {r["survey"] for r in back} == {"47B", "46B"}
    assert set(review.STATUS_COLUMNS) <= set(back[0])


def test_tool_written_fingerprints_collects_both_columns(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector" / "35_04_077").mkdir(parents=True)
    review.write_status("35_04_077", ROWS)
    fps = review.tool_written_fingerprints("35_04_077")
    assert "aa" * 32 in fps and "bb" * 32 in fps and "" not in fps


def test_worklist_is_one_file_for_the_whole_corridor(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    (tmp_path / "FMB_Vector").mkdir(parents=True)
    p = review.write_worklist({"35_04_077": ROWS, "35_04_074": []})
    with p.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2 and rows[0]["village_code"] == "35_04_077"


def _tracker_copy(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    return work


def test_the_auto_column_writer_is_retired():
    with pytest.raises(RuntimeError):
        review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "green"}])


def test_qc_passed_keeps_the_dropdowns_formatting_and_columns(tmp_path):
    import re
    work = _tracker_copy(tmp_path)
    before = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    result, missing = review.mark_qc_passed("35_04_074", ["39A"], workbook=work)
    assert result in ("applied", "no change") and missing == []
    after = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    assert after.count("<x14:dataValidation ") == before.count("<x14:dataValidation ")
    assert before.count("<cfRule") == after.count("<cfRule")
    head = lambda x: re.findall(r'<c r="([A-Z]+)1"', re.search(r'<row r="1".*?</row>', x, re.S).group(0))
    assert head(after) == head(before), "no column is ever added"


def test_qc_passed_changes_only_that_status_cell(tmp_path):
    import openpyxl
    work = _tracker_copy(tmp_path)
    rows = lambda: [list(r) for r in openpyxl.load_workbook(work)["Tracker"].iter_rows(values_only=True)]
    before = rows()
    target = [i for i, r in enumerate(before) if str(r[6]) == "74" and str(r[7]) == "39A"]
    assert target, "Thirukatchur 39A row exists"
    review.mark_qc_passed("35_04_074", ["39A"], workbook=work)
    after = rows()
    diffs = [(i, j) for i, (b, a) in enumerate(zip(before, after)) for j, (x, y) in enumerate(zip(b, a)) if x != y]
    assert all(i == target[0] and j == 9 for i, j in diffs)
    assert after[target[0]][9] == "QC Passed"


def test_qc_passed_reports_a_locked_workbook_and_unknown_surveys(tmp_path):
    work = _tracker_copy(tmp_path)
    assert review.mark_qc_passed("35_04_074", ["NOPE"], workbook=work) == ("no change", ["NOPE"])
    (tmp_path / "~$Tracker.xlsx").write_bytes(b"lock")
    assert review.mark_qc_passed("35_04_074", ["39A"], workbook=work)[0] == "locked"
