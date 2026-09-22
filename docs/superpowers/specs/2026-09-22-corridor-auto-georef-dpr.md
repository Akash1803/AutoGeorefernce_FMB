# Detailed Plan Report: auto-georeferencing the Puvi cadastral vectors of the rail corridor villages

Date 2026-09-22. Status: proposal for Akash's approval. No code has been changed for this plan.
Every number below was measured today on the files under `D:\Projects\Tambaram_Chengalpattu`,
`D:\Data\Puvi_data` and the repository at commit 8a27b3e, unless marked *estimate*.

## 1. The task as I read it

The Puvi archive gives, for each of the 21 revenue villages along the Tambaram - Chengalpattu
railway, a vector layer of whole survey numbers and a set of FMB sheets it has already
georeferenced. Both are only roughly on the ground. The task is to produce, for every survey in
scope, an FMB parcel placed the way you place one by hand: rigid (rotation and shift only, every
printed length kept), laid edge for edge against its neighbours, checked against satellite,
delivered as red outlines for your review and installed only on your word.

The Puvi vector is not the thing being moved. It is the survey list, a coarse locator, and the
fabric that says which parcels neighbour which. The product is a corrected parcel set that can
replace it.

Nothing is deleted or moved during this work. Every output is a new file beside the existing ones,
and the only files that change are the ones you say to install.

**One decision only you can take: what "whole vector data" means.**

| Scope | Surveys | Sheets in hand | Sheets to download |
|---|---|---|---|
| A. Everything touching the 30 m rail buffer (your 677-plot layer) | 677 | 657 | about 20 |
| B. Every survey of the 21 villages | 5,762 | 657 | about 5,100 |

This report plans scope A in full and shows what scope B adds. My recommendation is A first;
B is the same machinery run on more sheets, and its cost is dominated by downloading and
converting sheets, not by georeferencing.

## 2. What exists today

### 2.1 Data inventory

The corridor is 30.8 km of line. In rail order from the Tambaram end:

| Village | Code | km | Puvi surveys | In buffer | Sheets converted | Hand-placed | Automation so far |
|---|---|---|---|---|---|---|---|
| Peerkkararanai | 35_05_136 | 2.4-3.4 | 199 | 36 | 36 | 0 | - |
| Perunkalathur | 35_05_134 | 3.5-4.4 | 463 | 13 | 11 | 0 | - |
| Vandalur | 35_15_002 | 4.6-6.8 | 312 | 54 | 53 | 0 | - |
| Kelampakkam | 35_15_003 | 7.0-8.2 | 225 | 43 | 43 | 0 | - |
| Oorapakkam | 35_15_004 | 8.3-9.5 | 208 | 64 | 64 | 0 | - |
| Nanthivaram | 35_15_006 | 9.5-10.5 | 541 | 39 | 39 | 0 | - |
| Guduvancheri | 35_15_005 | 10.4-12.2 | 188 | 39 | 37 | 0 | - |
| **Thailavaram** | 35_04_052 | 12.3-13.5 | 94 | 44 | 52 | 55 (merged layer) | topology fixed today |
| **Potheri** | 35_04_054 | 13.5-14.8 | 243 | 22 | 20 | 0 | - |
| **Kattankulathur** | 35_04_056 | 15.1-17.7 | 362 | 38 | 38 | 0 | - |
| **Peramanur** | 35_04_076 | 17.9-19.2 | 77 | 8 | 9 | 0 | - |
| **Kizhikaranai** | 35_04_077 | 19.3-19.9 | 177 | 15 | 15 | 15 | ring runs, all 15 scored |
| **Thirukatchur** | 35_04_074 | 19.8-22.4 | 642 | 70 | 70 | 39 | 5 red awaiting your word, 26 left |
| Vinchiyambakkam | 35_04_085 | 22.6-23.0 | 31 | 10 | 10 | 0 | - |
| Thirutheri | 35_04_086 | 23.2-24.2 | 132 | 50 | 50 | 0 | - |
| Chettypunniyam | 35_04_071 | 24.2-25.6 | 462 | 36 | 36 | 0 | - |
| Paranur | 35_04_221 | 25.7-27.0 | 223 | 30 | 30 | 0 | - |
| Rajakulipettai | 35_04_224 | 27.0-28.3 | 68 | 39 | 38 | 0 | - |
| Pulipakkam | 35_04_225 | 28.3-29.1 | 287 | 6 | 6 | 0 | - |
| Natham A | 35_04_226 | 29.2-29.7 | 401 | 9 | 0 | 0 | - |
| Kolathur | 35_04_073 | not touching | 427 | 0 | 0 | 0 | not on the portal |

Bold rows are the stretch you named. Sheets are under `FMB_Sketches\<village>` as portal PDFs and
`FMB_Vector\<village>\<survey>_parcels.geojson` as converted vectors in sheet metres. Two-reader
neighbour transcription exists only for Kizhikaranai (15 sheets) and Thirukatchur (70 sheets).

Gaps between the buffer list and the sheets in hand: Potheri 226 and 54, Peramanur 43, Natham A's
nine, and a handful where the buffer count exceeds the sheet count (Perunkalathur, Vandalur,
Guduvancheri, Rajakulipettai). About 20 downloads in total. Survey naming differs between the
sources: Puvi writes `10 B`, the portal and the sheets write `10B`; 4 to 59 per village carry a
letter suffix, and up to 38 per village are written with a space. This needs one normalisation rule,
not hand work.

### 2.2 How good is Puvi as a starting position

Measured against your hand placements, survey by survey:

| Village | Matched | Centroid offset, median | 90th percentile | Max | Overlap with your parcel (IoU) median | Share with IoU over 0.5 |
|---|---|---|---|---|---|---|
| Thirukatchur, Puvi vector | 29 of 39 | 24.1 m | 52.4 m | 56.5 m | 0.27 | 14 % |
| Thirukatchur, Puvi georef FMB | 29 of 39 | 24.0 m | 52.4 m | 57.5 m | 0.25 | 10 % |
| Thailavaram, Puvi vector | 30 of 55 | 6.5 m | 11.1 m | 17.7 m | 0.73 | 97 % |
| Thailavaram, Puvi georef FMB | 27 of 55 | 6.6 m | 11.3 m | 56.3 m | 0.69 | 81 % |

Two conclusions. Puvi's own georeferenced FMB is no better than its vector, so there is nothing to
gain from it. And Puvi's error is a village-level property: 6 m in Thailavaram, 24 m in
Thirukatchur, with the shapes themselves right (area ratio median 1.00 in both). A parcel can be
found from Puvi; it cannot be placed from Puvi. This confirms the reference-only rule from the
Kizhikaranai work.

### 2.3 What the automation can do today

The engine (`autogeoref`, SOP version 3) places a sheet by laying it against neighbours already on
the ground, matching shared boundary runs edge for edge, ranking candidate poses by printed
neighbour numbers and sides, adjusting the block, and checking the result against Google
Satellite. It needs, per village: converted sheets, at least one hand-placed seed with its
`.points` file, and the two-reader transcription of the printed neighbour numbers.

Evidence from the harness (`_logs\eval\rows.csv`, 158 rows, 58 runs, Kizhikaranai and
Thirukatchur):

- 48 distinct parcels have been placed by the tool: 1 green, 29 amber, 16 red.
- 20 parcel-runs were scored against a hand placement: median 4.1 m, 90th percentile 7.1 m,
  worst 15.5 m. Only 14 distinct labelled parcels exist, so this is indicative, not statistical.
- Your own placements disagree with each other by 1 to 3 m (anchor fit rms 0.19 to 3.05 m).
  That is the accuracy floor: nothing can be scored finer than that.
- Boundary fit between placed neighbours: median 0.5 m rms.
- A 15-parcel village takes about 10 minutes of machine time; the first run per village downloads
  about 65 MB of satellite.

Known limits that matter for a corridor-scale run:

- Strip parcels (railway land, roads, tanks) cannot be placed from their end edges and are hand
  seeds by rule.
- Chain matching is position-blind: a wrong placement can score perfectly inside the tool. Only
  hand placements and satellite catch it. This is the reason for the Puvi window in section 4.
- Sheets can be wrong: 613 in Thirukatchur was drawn 9.9 m short and broke four neighbours until
  it was corrected. Sheet QC exists (`--sheet-qc`) and was run today on the six stretch villages:
  Thailavaram 4 sheets flagged, Potheri 4, Kattankulathur 4, Peramanur 2 of 9 (43A and 43B are
  long thin water or road strips), Kizhikaranai 3, Thirukatchur 4 (526A alone has 177 subdivision
  plots without numbers).
- The engine works one village at a time and reads its village code from the folder. Parcels on a
  village boundary have neighbours in the next village that it cannot see.
- **Transcription is the binding constraint.** Without the two-reader neighbour table the printed
  contradiction and side checks are silent, and the tool once matched 569C to 609, 1.1 km away, in a
  village where only four sheets had been read. Seventeen of nineteen villages have no table yet.
- A Puvi-seeded pose search (`window.start_poses`: 2° then 0.5° rotation steps, 2 m shift steps
  over a 6 m range) already exists in the code with tests but no caller. Its 6 m range was written
  for Kizhikaranai and is a quarter of Thirukatchur's Puvi error; the seedless start in section 4
  builds on it with the calibrated window.
- Scaling walls for one big run, read from the code: every sheet is re-read from disk on every
  use (no cache); growth is passes over all pending against all placed, cubic in the worst case;
  the satellite export holds the whole window in memory (a 3 km window is 1.2 GB) and refuses the
  export if over 2 % of tiles come back empty; every output file makes a round trip to QGIS to see
  whether the layer is held; the review group loads one QGIS layer per parcel. Measured: a ring-2
  run took 928 s, leave-one-out about 40 s per parcel, Thirukatchur's automatic window with 19
  anchors ran past ten minutes.
- Compute is your laptop only: no GPU, no cloud, Google tiles for inference only, nothing trained.

## 3. Options

### Option A: the whole corridor in one run

Run all 19 villages with sheets as one job, using Puvi as the seed for every parcel.

Why not: Puvi is 6 to 57 m off, the tool cannot tell a wrong placement from a right one without a
hand-placed neighbour or satellite, no transcription exists for 17 villages, and the code as it
stands would not finish: the imagery export alone exceeds memory at corridor size and the growth
passes are cubic. One run would
produce 657 red parcels of unknown quality and a review load nobody can clear. It also hides the
answer to the only question that matters, which is whether the tool can start a village that has no
hand seed.

### Option B: pilot on the stretch, then roll out village by village (recommended)

Use the Thailavaram - Thirukatchur stretch, because it is the only place where both ends are
already on the ground by your hand: Thailavaram (55 surveys) at km 12.3 and Kizhikaranai plus
Thirukatchur (54 surveys) at km 19.3. Potheri, Kattankulathur and Peramanur sit between them with
no seed at all, which is exactly the situation the other 13 villages are in. The pilot answers, on
68 sheets, whether a village can be started from a Puvi-located, satellite-confirmed seed instead
of a hand-placed one, and it finishes Thirukatchur's 26 remaining sheets on the way.

The pilot's gate decides the roll-out. If it passes, the remaining 13 villages run in rail order in
batches of two or three, each batch reviewed by you before the next.

### Option C: village by village with a hand seed each, as now

Continue the current SOP: you place one or two seeds per village by hand, the tool grows rings.

This works and is the fallback if option B's seedless start fails. Its cost is your time: 17
villages times two seeds in the Georeferencer, plus two-reader transcription of each village. It
does not scale to scope B (5,762 surveys).

## 4. Recommended plan

### Phase 0: preparation for every village (no georeferencing yet)

1. **Survey-name reconciliation.** One rule maps Puvi `10 B`, portal `10B` and sheet `10B` to one
   key; output `survey_index_<village>.csv` with Puvi polygon, sheet file, in-buffer flag, and
   `missing sheet` where there is none. Kolathur is dropped with a note; Natham A is queued for
   download.
2. **Missing sheets.** Download the roughly 20 in-buffer surveys with no sheet through the FMB
   downloader, convert them with the existing batch, enrich attributes. Known portal rules apply:
   numeric numbers first, letter suffixes as units when the plain number fails, strip sheets at
   scale 2000.
3. **Sheet QC on all 19 villages** (done for six today). Each flagged sheet is looked at once;
   a drawing error is recorded in `sheet_corrections.csv` as for 613, never edited silently.
4. **Puvi calibration per village.** Where hand placements exist, the offset table of section 2.2
   is the calibration. Where none exist, the village inherits the corridor-wide figure (median 24 m,
   90th percentile 52 m) until its first accepted parcels tighten it. This number becomes the
   search window for the seedless start.
5. **Cross-village fabric.** Build the neighbour graph over the whole corridor, not per village, so
   a parcel on the Thailavaram - Potheri boundary sees its neighbour across the line. The engine
   already takes any placed body as a settled neighbour; the change is in what is loaded, not how it
   matches.

Deliverable: one corridor index, sheets complete, QC list, calibration table. Nothing on the map
changes.

### Phase 1: pilot, Thailavaram to Thirukatchur

Sheets in play: Potheri 20 (+2 to download), Kattankulathur 38, Peramanur 9, Thirukatchur's 26
unplaced. Fixed truth: Thailavaram 55, Kizhikaranai 15, Thirukatchur 39, none of them ever moved.

1. **Transcription** of the printed neighbour numbers and sides for the 68 new sheets, two
   independent readers as before, merged with the review CSV you resolve. Thirukatchur's 70 are
   already read.
2. **Seedless start, the new piece.** For a village with no hand seed, the tool takes the Puvi
   polygon of each candidate parcel as a prior pose with the calibrated uncertainty, generates
   candidate poses inside that window from the sheet outline against satellite edges, and admits a
   seed only when three conditions hold: satellite share above the green threshold, a runner-up
   pose that explains under 80 % of the best, and printed neighbours that agree with the Puvi
   fabric on every side. A parcel that fails is not seeded; the village waits for the next best
   candidate or for a hand seed from you. Expected: two or three seeds per village.
3. **Ring growth** from each accepted seed as in the SOP, rings reviewed by you, accepted parcels
   becoming fixed. Cross-village: Potheri's first ring can grow from Thailavaram's hand-placed
   parcels at km 13.5 without any seed of its own; Peramanur's from Kizhikaranai.
4. **Held-out scoring.** Before Thirukatchur's 26 are placed, the tool re-places 13 of your 39
   hand parcels with those 13 withheld (leave-one-out, as on 2026-09-18) so the pilot has a
   number that is not self-reported.
5. **Topology** on the automated parcels only, with the guard rules from this week (carved
   parcel goes red, no fills beyond 25 m², no needles).
6. **Delivery**: red outlines in your project at your line width, `georef_status.csv`, GCP
   `.points` per parcel, tracker columns only if you ask.

**Gate to pass before phase 2:**

| Measure | Pass |
|---|---|
| Held-out Thirukatchur parcels within 3 m of your placement | at least 10 of 13 |
| Seedless start in Potheri and Kattankulathur | at least one accepted seed per village, confirmed by you |
| Parcels you reject at review | under 20 % of those delivered |
| Wrong placements that came out green | zero |

If the seedless start fails but ring growth from the neighbouring village works, the roll-out
proceeds with cross-village growth and hand seeds only where a village has no placed neighbour
(option C for those few).

### Phase 2: roll-out, remaining 13 villages, 430 sheets

In rail order, batches of two or three villages, each batch: transcription, seedless start or
growth from the previous batch's accepted parcels, rings, your review, install on your word.
Guduvancheri - Nanthivaram - Oorapakkam first (they touch Thailavaram's end of the stretch), then
westward to Peerkkararanai, then east from Vinchiyambakkam to Natham A.

### Phase 3, only if scope B is chosen

The same pipeline on the 5,100 surveys outside the buffer. New cost is almost entirely downloads,
conversion and transcription. Georeferencing inside a village becomes easier, not harder, as the
placed fabric fills in. Not planned in detail until scope A is done.

## 5. Engineering work the plan needs

Each item is one PR, reviewed before the next, in the order the phases need them.

| # | Change | Where | Why |
|---|---|---|---|
| 1 | Survey-key normalisation and corridor index | new `autogeoref/corridor.py`, `neighbours.normalise` | one key across Puvi, portal, sheets |
| 2 | Puvi prior: per-village offset calibration, search window per parcel | `references.py`, `config.py` | the seedless start needs a bounded window |
| 3 | Seedless seed selection with the three admission conditions | `engine.py` (new `seed_candidates`), `align.py` | start a village without a hand seed |
| 4 | Cross-village settled bodies: load neighbours from adjacent villages as anchors | `anchors.py`, `paths.py` | boundary parcels |
| 5 | Batch runner over villages with per-batch logs and a corridor status table | `georef_village.py` or new `georef_corridor.py` | 19 villages without hand-driving each |
| 6 | Harness: corridor-level report, held-out scoring as a command | `evalrows.py`, `report.py` | the gate needs numbers |
| 7 | Throughput: cache parsed sheets in memory, expose ring runs, `only` and `use_imagery` on the CLI, one QGIS round trip per run instead of per file, review group per batch | `sheets.py`, `georef_village.py`, `files.py`, `review.py` | a 40-sheet village must run in minutes, not hours |
| 8 | Imagery per ring, not per village: keep the local 300 m window rule of the SOP and never export a whole village | `engine.py`, `raster.py` | memory wall and the 2 % empty-tile refusal |

No new dependency. No training. Google tiles remain inference-only and the cache stays deletable.

## 6. Effort and time (estimates, marked as such)

| Item | Basis | Estimate |
|---|---|---|
| Phase 0 | 19 villages, scripts exist for QC and conversion | 2 working days |
| PRs 1-3 | code plus tests, one at a time with your review | 3 to 4 working days |
| Pilot transcription, 68 sheets, two readers | Thirukatchur's 70 took one session | 1 day |
| Pilot runs and your reviews | ring runs of about 15 min each, held-out scoring about 40 s per parcel; review is the long pole | 3 to 5 days elapsed |
| PRs 4-8 | items 7 and 8 are plumbing, not algorithms | 3 to 4 days |
| Phase 2, 430 sheets in 5 or 6 batches | transcription and your review dominate | 3 to 4 weeks elapsed |

Machine time is negligible. The schedule is set by transcription and by how fast you can review.

## 7. Risks

| Risk | Effect | Mitigation |
|---|---|---|
| Seedless start places a parcel confidently in the wrong place (position-blind matching) | wrong fabric grows from it | three admission conditions; Puvi window caps the error at the village calibration; your review of every seed before any ring |
| Puvi error larger than calibrated in an unmeasured village | seed search misses | window widens to the corridor 90th percentile (52 m) on a miss, never silently |
| Sheet drawn wrong (as 613) | neighbours break | sheet QC first; corrections recorded, never silent edits |
| Strip parcels (rail, road, tank) | unplaceable by edges | hand seeds by rule, flagged in the corridor index up front |
| Name mismatch between sources | parcel matched to the wrong survey | normalisation PR with tests on the odd forms seen (`1 A`, `10 B`, `569/C`) |
| Review load | the plan stalls | batches of 2-3 villages; nothing installs without your word, so nothing waits on the machine |
| Scope B downloads (5,100 sheets) | weeks of portal time | not started until scope A is accepted |

## 8. What you get at the end of scope A

- For each of the 677 in-buffer surveys: a placed parcel file, its GCPs, its status row, and a
  colour; installed as your `_parcels_modified.gpkg` only where you said so.
- A corridor status table and an accuracy report against every hand placement that exists.
- Puvi untouched, every original file kept, every run logged under `_logs`.

## 9. Questions for you

1. Scope A (buffer, 677) or scope B (whole villages, 5,762)?
2. Is the pilot stretch as I read it right: Potheri, Kattankulathur and Peramanur are the seedless
   test, Thirukatchur's 26 get finished, Thailavaram and Kizhikaranai stay as they are?
3. The gate numbers in phase 1: are 3 m and 10 of 13 the bar you want, or should it be tighter?
