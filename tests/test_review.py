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


def test_tracker_update_keeps_the_dropdowns_and_formatting(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    before = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "green",
                            "notes": "anchored", "run": "2026-09-18T10:00:00"}], workbook=work)
    after = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    assert after.count("<x14:dataValidation ") == before.count("<x14:dataValidation ")
    assert "<conditionalFormatting" in after
    assert before.count("<cfRule") == after.count("<cfRule")


def test_tracker_update_never_writes_the_status_column(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    import openpyxl
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    before = [r[9] for r in openpyxl.load_workbook(work, data_only=True)["Tracker"].iter_rows(values_only=True)]
    review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "red",
                            "notes": "x", "run": "2026-09-18T10:00:00"}], workbook=work)
    after = [r[9] for r in openpyxl.load_workbook(work, data_only=True)["Tracker"].iter_rows(values_only=True)]
    assert before == after, "Status is the team's column"


def test_tracker_update_reports_a_locked_workbook_instead_of_failing(tmp_path):
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    (tmp_path / "~$Tracker.xlsx").write_bytes(b"lock")
    assert review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "green",
                                   "notes": "x", "run": "r"}], workbook=work) == "locked"


def test_tracker_filter_reaches_the_new_columns(tmp_path):
    """526A was invisible on 2026-09-18 because the autofilter ref stopped short of its row."""
    src = paths.TRACKER
    if not src.exists():
        pytest.skip("tracker not present")
    work = tmp_path / "Tracker.xlsx"
    work.write_bytes(src.read_bytes())
    review.update_tracker([{"village_code": "35_04_077", "survey": "47B", "colour": "green",
                            "notes": "x", "run": "r"}], workbook=work)
    sheet = zipfile.ZipFile(work).read("xl/worksheets/sheet1.xml").decode("utf-8")
    import re
    flt = re.search(r'<autoFilter ref="A1:([A-Z]+)\d+"', sheet).group(1)
    dim = re.search(r'<dimension ref="A1:([A-Z]+)\d+"/>', sheet).group(1)
    assert flt == dim, "the filter must cover every column, not just the team's"
    assert review._index(flt) == review._index("O") + len(review.TRACKER_COLUMNS)
