"""The evaluation report: comparison table, worst ten, and the sentence that goes on every one.

Metrics are computed per village fold and never pooled, because the split is by village. Each
fold is reported twice: all labelled rows, and without the middle stretch tier, so the effect of
keeping that tier in the labels is visible.
"""
import datetime
import logging
import statistics
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from . import evalrows

log = logging.getLogger(__name__)

BOILERPLATE = ("Errors are measured against the team's hand placements, which disagree with each "
               "other by 1 to 3 m. Results are therefore capped by anchor consistency and are not "
               "absolute accuracy. Stretched parcels are scored against both the FMB-exact and the "
               "saved placement and carry no label.")

CAUSES = [("contradict", "printed side or neighbour contradicted"),
          ("does not reach printed", "does not reach a printed neighbour"),
          ("no placed neighbour", "no placed neighbour"),
          ("outside the satellite window", "outside the satellite window"),
          ("next best", "runner-up pose close"),
          ("railway", "clipped by railway land"),
          ("backfilled", "backfilled row, features missing")]


def _f(x: object) -> Optional[float]:
    try:
        return None if x in (None, "") else float(x)
    except (TypeError, ValueError):
        return None


def _pct(num: int, den: int) -> Optional[float]:
    return round(100.0 * num / den, 1) if den else None


def _q(values: List[float], q: float) -> Optional[float]:
    if not values:
        return None
    vs = sorted(values)
    k = (len(vs) - 1) * q
    lo, hi = int(k), min(int(k) + 1, len(vs) - 1)
    return round(vs[lo] + (vs[hi] - vs[lo]) * (k - lo), 2)


def metrics(rows: Iterable[Dict[str, object]], acceptance_m: float = 3.0) -> Dict[str, object]:
    """Counts and rates for one set of rows; rule-based colour against the label.

    Sheet-QC failure rows carry a reason and no label and are left out of every figure here.
    """
    rows = [r for r in rows if not evalrows.is_qc_failure(r)]
    lab = [r for r in rows if evalrows.is_labelled(r)]
    within = [r for r in lab if str(r.get("within_3m")) == "1"]
    green = [r for r in lab if r.get("colour") == "green"]
    red = [r for r in lab if r.get("colour") == "red"]
    green_amber_within = [r for r in within if r.get("colour") in ("green", "amber")]
    err_fmb = [v for v in (_f(r.get("err_fmb_m")) for r in rows) if v is not None]
    err_hand = [v for v in (_f(r.get("err_hand_m")) for r in rows) if v is not None]
    return {
        "rows": len(rows), "labelled": len(lab), "within_pct": _pct(len(within), len(lab)),
        "green": len(green), "green_within_pct": _pct(sum(1 for r in green if str(r.get("within_3m")) == "1"), len(green)),
        "red": len(red), "red_beyond_pct": _pct(sum(1 for r in red if str(r.get("within_3m")) == "0"), len(red)),
        "recall_green_amber_pct": _pct(len(green_amber_within), len(within)),
        "err_fmb_median_m": _q(err_fmb, 0.5), "err_fmb_p90_m": _q(err_fmb, 0.9),
        "err_hand_median_m": _q(err_hand, 0.5), "err_hand_p90_m": _q(err_hand, 0.9),
        **{k: v for k, v in evalrows.distinct_summary(rows).items() if k.startswith(("distinct", "parcels"))},
    }


PR3_MIN_PARCELS = 300
PR3_MIN_VILLAGES = 2


def pr3_gate(rows: List[Dict[str, object]], min_parcels: int = PR3_MIN_PARCELS, min_villages: int = PR3_MIN_VILLAGES) -> Dict[str, object]:
    """The PR 3 gate counts distinct labelled parcels across villages, never rows."""
    current = [r for r in rows if str(r.get("backfilled", "0")) in ("", "0")]
    d = evalrows.distinct_summary(current)
    ok = d["distinct_labelled"] >= min_parcels and len(d["villages"]) >= min_villages
    return dict(d, ok=ok, min_parcels=min_parcels, min_villages=min_villages)


def comparison_table(rows: List[Dict[str, object]], acceptance_m: float = 3.0) -> List[Dict[str, object]]:
    """One line per village, per row set (all labelled / without the middle tier), per method."""
    out: List[Dict[str, object]] = []
    villages = sorted({str(r.get("village")) for r in rows})
    for v in villages:
        vrows = [r for r in rows if str(r.get("village")) == v]
        current = [r for r in vrows if str(r.get("backfilled", "0")) in ("", "0")]
        band_free = [r for r in current if not str(r.get("stretch_band", "")).startswith("5-")]
        # backfilled rows come from runs of older code and are kept for the record, not as the baseline
        for set_name, subset in (("current code, all labelled", current),
                                 ("current code, without 5-10 % tier", band_free),
                                 ("all rows incl. backfilled", vrows)):
            m = metrics(subset, acceptance_m)
            out.append(dict(village=v, row_set=set_name, method="rule-based", **m))
            if any(r.get("p_within") not in (None, "") for r in subset):
                out.append(dict(village=v, row_set=set_name, method="candidate", **metrics(subset, acceptance_m)))
    return out


def cause_from_notes(notes: str) -> str:
    text = (notes or "").lower()
    for key, label in CAUSES:
        if key in text:
            return label
    return "see notes"


def worst10(rows: List[Dict[str, object]], include_backfilled: bool = False) -> List[Dict[str, object]]:
    """The ten largest errors against the FMB-exact reference, current code only unless asked."""
    if not include_backfilled:
        rows = [r for r in rows if str(r.get("backfilled", "0")) in ("", "0")]
    rows = [r for r in rows if not evalrows.is_qc_failure(r)]
    scored = [(r, _f(r.get("err_fmb_m"))) for r in rows]
    scored = [(r, e) for r, e in scored if e is not None]
    scored.sort(key=lambda re: -re[1])
    out = []
    for r, e in scored[:10]:
        out.append({"village": r.get("village"), "survey": r.get("survey"), "run_id": r.get("run_id"),
                    "err_fmb_m": e, "err_hand_m": _f(r.get("err_hand_m")), "colour": r.get("colour"),
                    "stretched": r.get("stretched"), "cause": cause_from_notes(str(r.get("notes") or ""))})
    return out


def _md_table(lines: List[Dict[str, object]], columns: List[str]) -> str:
    head = "| " + " | ".join(columns) + " |\n|" + "|".join("---" for _ in columns) + "|\n"
    body = "".join("| " + " | ".join("" if l.get(c) is None else str(l.get(c)) for c in columns) + " |\n" for l in lines)
    return head + body


def write_report(rows: List[Dict[str, object]], out_path: Path, run_label: str,
                 acceptance_m: float = 3.0, extra: Optional[str] = None) -> Path:
    """Markdown report with the boilerplate first, then the table, then the worst ten."""
    table = comparison_table(rows, acceptance_m)
    worst = worst10(rows)
    cols = ["village", "row_set", "method", "rows", "labelled", "distinct_parcels", "distinct_labelled",
            "parcels_within", "parcels_beyond", "within_pct", "green", "green_within_pct",
            "recall_green_amber_pct", "red", "red_beyond_pct", "err_fmb_median_m", "err_fmb_p90_m",
            "err_hand_median_m", "err_hand_p90_m"]
    gate = pr3_gate(rows)
    qc_rows = [r for r in rows if evalrows.is_qc_failure(r)]
    wcols = ["village", "survey", "run_id", "err_fmb_m", "err_hand_m", "colour", "stretched", "cause"]
    text = ["# Evaluation report: %s" % run_label,
            "", datetime.datetime.now().isoformat(timespec="seconds"), "",
            BOILERPLATE, "",
            "Acceptance: within %.1f m of the FMB-exact reference. Split: by village; no pooled figure." % acceptance_m,
            "", "## Comparison", "", _md_table(table, cols),
            "Percentages: `within_pct` labelled rows within acceptance; `green_within_pct` share of green rows within; "
            "`recall_green_amber_pct` share of within-acceptance rows ranked green or amber; `red_beyond_pct` share of red rows beyond.",
            "", "## Labelled parcels and the PR 3 gate", "",
            "Distinct labelled parcels (latest row per parcel): %d, of which %d within and %d beyond acceptance, "
            "across %d village(s): %s. Rows are not the gate; parcels are. Gate (%d parcels, %d villages): %s."
            % (gate["distinct_labelled"], gate["parcels_within"], gate["parcels_beyond"], len(gate["villages"]),
               ", ".join(gate["villages"]) or "-", gate["min_parcels"], gate["min_villages"], "MET" if gate["ok"] else "not met"),
            "", "Sheet-QC failures logged with a reason and no label, excluded from every figure above: %d (%s)."
            % (len(qc_rows), ", ".join("%s %s: %s" % (r.get("village"), r.get("survey"), r.get("sheet_qc_reason")) for r in qc_rows) or "none"),
            "", "## Worst 10 by error against the FMB-exact reference", "", _md_table(worst, wcols)]
    if extra:
        text += ["", extra]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(text) + "\n", encoding="utf-8")
    log.info("report written to %s", out_path)
    return out_path
