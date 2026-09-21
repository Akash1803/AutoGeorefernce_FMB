from autogeoref import report


def _row(village, survey, colour, err, band="<=5%", label=None, notes="", err_hand=None):
    if label is None:
        label = "" if err is None else ("1" if err <= 3.0 else "0")
    return {"village": village, "survey": survey, "run_id": "r1", "colour": colour, "err_fmb_m": err,
            "err_hand_m": err if err_hand is None else err_hand, "stretch_band": band, "within_3m": label,
            "notes": notes, "stretched": "0"}


ROWS = [_row("A", "1", "green", 1.0), _row("A", "2", "green", 4.0), _row("A", "3", "amber", 2.0),
        _row("A", "4", "red", 16.0), _row("A", "5", "amber", 2.5, band="5-10%"),
        _row("A", "6", "red", 2.0), _row("B", "1", "amber", 1.0),
        _row("B", "2", "amber", 9.0, band="5-10%", label="", notes="does not reach printed neighbour(s) 42B")]


def test_metrics_count_what_the_colour_rule_gets_right():
    m = report.metrics([r for r in ROWS if r["village"] == "A"])
    assert m["rows"] == 6 and m["labelled"] == 6
    assert m["green"] == 2 and m["green_within_pct"] == 50.0
    assert m["red"] == 2 and m["red_beyond_pct"] == 50.0
    assert m["recall_green_amber_pct"] == 75.0            # within: 1,3,5,6 -> green/amber: 1,3,5
    assert m["err_fmb_median_m"] == 2.25 and m["err_fmb_p90_m"] > 4.0


def test_comparison_table_reports_each_village_with_and_without_the_middle_tier():
    t = report.comparison_table(ROWS)
    assert [(l["village"], l["row_set"]) for l in t] == [("A", "all labelled"), ("A", "without 5-10 % tier"),
                                                          ("B", "all labelled"), ("B", "without 5-10 % tier")]
    assert t[0]["labelled"] == 6 and t[1]["labelled"] == 5
    assert all(l["method"] == "rule-based" for l in t)


def test_worst10_is_sorted_and_names_a_cause():
    w = report.worst10(ROWS)
    assert [(x["village"], x["survey"]) for x in w[:3]] == [("A", "4"), ("B", "2"), ("A", "2")]
    causes = {(x["village"], x["survey"]): x["cause"] for x in w}
    assert causes[("B", "2")] == "does not reach a printed neighbour"
    assert causes[("A", "2")] == "see notes"


def test_report_file_starts_with_the_accuracy_cap_sentence(tmp_path):
    p = report.write_report(ROWS, tmp_path / "r.md", "unit test")
    text = p.read_text(encoding="utf-8")
    assert report.BOILERPLATE in text
    assert text.index(report.BOILERPLATE) < text.index("## Comparison")
    assert "Worst 10" in text and "| A | all labelled |" in text
