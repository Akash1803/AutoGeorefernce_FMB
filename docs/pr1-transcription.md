# PR 1: two-reader neighbour transcription

Status: implemented 2026-09-21 (`autogeoref/transcribe.py`), first run on Thirukatchur 35_04_074.
Follows the approved [ml-hybrid-plan.md](ml-hybrid-plan.md) and Akash's ordered steps of
2026-09-21: PR 0 merged first; run on the village; report two-reader agreement and write every
disagreement to a review CSV; he checks the disagreements plus ten random agreed sheets against
the PDFs; accuracy per field from that sample; then the minimum seed set from the neighbour
graph; seeds are anchors and are never scored.

## What it does

1. **Render.** `georef_village <village> --render-sheets` crops the drawing area of every sheet
   PDF to `<project>/_cache/sheet_renders/<village>/<survey>.png` (gitignored, deletable; scale
   `transcription.render_scale`, default 1.6 px per PDF point).
2. **Read.** Two readers each look at every crop and write
   `<vector_dir>/neighbour_transcription/readings/<name>.json`:
   `{"<survey>": [{"number": "569/C", "side": "NE", "confidence": "high"}, ...]}`.
   A label is a black underlined survey number outside the outline; `side` is the compass side
   of the label around the outline; blue plot numbers, red corner letters, dimensions and
   `R.S. No.` stones are not labels. The backend names come from `transcription.readers`
   (`json:A`, `json:B`); a reader is any vision model session, in a separate context per reader.
   No OCR runs in this PR; PaddleOCR (PR 4) will be another backend name.
3. **Merge.** `georef_village <village> --transcribe`:
   - numbers are compared normalised (`569/C` = `569C`, village labels dropped);
   - a number both readers print is trusted (`readers: 2`); its side is the one both gave, or
     empty with a `side` disagreement;
   - a number one reader prints stays untrusted (`readers: 1`, `confidence: low`, the engine
     ignores it) and becomes a `number` disagreement;
   - writes `nb_override_<village>.json` (the engine's table), `review_<village>.csv` (one row per
     disagreement, `resolution` column for the analyst: yes/no for a number row, a compass side
     for a side row; resolutions survive re-runs), `agreement_<village>.md` (agreement per field,
     the disagreements, ten random fully-agreed sheets to check, per-sheet table) and
     `seeds_<village>.json`.
4. **Seeds.** From the trusted graph (edge when either sheet prints the other, restricted to
   sheets we have), a greedy dominating set: every sheet is a seed or touches one. Two answers:
   from scratch (ties broken towards parcels already placed) and additional to what is placed.
   Labelled parcels = placed parcels that are not seeds.

## Accuracy

Both readers used on 2026-09-21 are the same vision model in separate contexts. Their agreement
therefore overstates accuracy against the PDF. The accuracy figure per field (number, side) comes
from Akash's check of every disagreement plus the ten sampled agreed sheets; it is reported
after that check, not before.

## Files

| File | Purpose |
|---|---|
| `autogeoref/transcribe.py` | render, readers, merge, review CSV, report, graph, seeds |
| `configs/default.json` `transcription` | readers, render dir and scale, sample size and seed |
| `tests/test_transcribe.py` | merge rules, resolutions round trip, dominating set, render, config |
| `<vector_dir>/neighbour_transcription/readings/*.json` | the readers' raw output (project data) |
| `<vector_dir>/neighbour_transcription/review_<village>.csv` | disagreements and resolutions |
| `<vector_dir>/neighbour_transcription/agreement_<village>.md` | the agreement report |
| `<vector_dir>/neighbour_transcription/seeds_<village>.json` | the seed plan |

Every engine run now logs evaluation rows by default (`--no-rows` to skip).

## Follow-ups of 2026-09-21 (Akash's review of PR 0 and PR 1)

- Every table entry carries `source`: `reader` (both readers agreed) or `hand` (set through the
  review CSV; such rows also carry `resolved_by = hand`). Reader agreement is computed before
  resolutions, so hand-resolved rows never count towards reader accuracy. The review CSV has a
  `resolved_by` column.
- Sheets whose outline was corrected by hand (`sheet_corrections.csv`, rows other than plot
  re-conversions, with `verified_by`) are excluded from the seed set and from labels; the seed
  plan lists them under `excluded_from_seeds_and_labels`. Today: 613 (`verified_by = eye`).
- The two readers on 2026-09-21 were the same vision model (Claude Fable 5.1) in two separate
  contexts with different sheet batches. Not two different models.
- Kizhikaranai (35_04_077) was re-run through this flow: 15 sheets, number agreement 1.0, side
  agreement 0.945, 4 side rows to review; trusted numbers identical to the earlier table.
  There is no Kolathur (35_04_073) in this project; the earlier data plan was wrong about it.
