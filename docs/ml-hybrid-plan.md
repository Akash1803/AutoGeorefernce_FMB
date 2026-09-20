# Hybrid ML/DL plan for FMB auto-georeferencing: Step 1 audit and Step 2 architecture

Date 2026-09-20. Status: proposal, no code changed. Written against commit 15147da of
`AutoGeorefernce_FMB`. Everything below is either read from the repository or measured on the
Kizhikaranai runs; nothing is invented. Where I could not verify a number I say so.

## 0. Current state, filled in

| Item | Value |
|---|---|
| Repo | `D:\code\FMB_to_GeoJSON` (GitHub Akash1803/AutoGeorefernce_FMB), 4071 lines in `autogeoref/` + `fmb_to_geojson.py` |
| Runtime | Python 3.12.10, OpenCV 5.0, rasterio 1.5, geopandas 1.1.4, shapely 2.1.2, scipy 1.18, numpy 2.5.1, scikit-learn 1.9. **No PyTorch installed, no GPU on this machine.** |
| Input | Tamil Nadu FMB sheets as vector PDFs from the Collabland portal (656 in the 60 m rail corridor, 19 villages), 1:500 to 1:2000 |
| Reference | The team's hand placements (32 GeoPackages: 15 Kizhikaranai, 11 Kolathur, 6 superseded; 26 `.points` files); Google Satellite XYZ tiles at zoom 20 (0.15 m/px) exported per window; Puvi village vectors (reference only, 22 identity conflicts known); the two-reader transcription of printed neighbour numbers and sides (Kizhikaranai only so far) |
| Target CRS | EPSG:32644 |
| Volume | 656 sheets once for the corridor, then per new corridor; no daily flow |
| Accuracy target | Not set by you yet. Measured floor: your own placements disagree with each other by 1 to 3 m (anchor fit rms 0.19 to 3.05 m). See question 1. |

## 1. Step 1: audit of the pipeline as it is

### 1.1 Stages, inputs and outputs

| # | Stage | Module | Input | Output | Method today |
|---|---|---|---|---|---|
| A | Sheet to vector | `fmb_to_geojson.py` | portal PDF (vector drawing + glyph text) | `<survey>_parcels.geojson` in sheet metres, per-plot attributes, review flags | PDF path extraction; glyph template matching against `glyphlib.json` (213 templates); size-group heuristics scaled off A0 height |
| B | Neighbour transcription | `neighbours.py` | printed numbers and sides on the sheet | `nb_override_<village>.json` | done by a two-reader agent workflow, not by code; code only reads and normalises it |
| C | Anchors | `anchors.py` | team GeoPackages + `.points` | pose (rotation, shift) per hand-placed parcel, fingerprint, rms, scale | rigid Procrustes fit of the `.points`; fallback ICP-style vertex fit on the geometry |
| D | Candidate poses | `engine.pose_candidates`, `match.common_chains` | free sheet outline + placed neighbour outlines | list of (theta, t) | chainage walk over edge lengths and turn angles, rigid fit per matched run |
| E | Candidate ranking | `engine.rank_poses`, `side_checks`, `printed_contradictions`, `_overlap` | candidates, placed bodies, transcription | ranked list, support, contradictions | hand-written score: contradictions first, then supported boundary length |
| F | Image evidence | `raster.py`, `edges.py`, `align.py` | Google tiles, top-3 poses | share, margin, observability, possibly an image pose | OpenCV LSD line segments, KD-tree, binary share metric, grid search over shift and rotation |
| G | Block adjustment | `fit.block_adjust` | pair, line, prior observations | adjusted poses, sigmas | scipy least squares, Huber loss, centred frame, 3 rounds |
| H | GCPs | `gcp.py` | corners + segment index | `.points` (QGIS Georeferencer format) + CSV | corner within 1.5 m of a segment = matched |
| I | Topology | `topology.py` | placed bodies, team bodies, rail set | clipped/filled plots, `rail_conflicts.csv` | clip overlaps, fill slivers, clean |
| J | Verdict | `engine.colour_of` | row features | green/amber/red | hand-calibrated rule on the 15-parcel leave-one-out |
| K | Review | `review.py` | rows | status CSV, worklist, QGIS group, tracker columns | XML patching, QGIS socket |

### 1.2 Every hard-coded threshold and assumption

Read from the source with `grep` on 2026-09-20; the comment column is what the code says the number is for.

| Module | Name | Value | Why it is that value |
|---|---|---|---|
| align | MATCH_DIST | 1.5 m | spike tolerance for a sample to count as on an image edge |
| align | MATCH_DEG | 12° | edge direction agreement |
| align | EXCL_M / EXCL_DEG | 3.0 m / 3.0° | exclusion zone when looking for a competing peak |
| align | MARGIN_MIN | 0.15 | best pose must beat the best outside pose by this share |
| align | PEAK_RATIO | 0.75 | local maxima above this fraction count as competitors |
| align | PRIMARY_MIN / CROSS_MIN / CROSS_ANGLE | 40 m / 20 m / 30° | observability: matched length in two directions |
| align | GROUP_TOL | 10° | bearings within this are one direction |
| anchors | POSE_RMS_LIMIT | 3.0 m | above this the `.points` fit is distrusted, geometry fit used |
| anchors | VERTEX_REACH | 4.0 m | how far a sheet vertex may look for its hand-edited twin |
| anchors | CONFLICT_SAMPLE | 1.0 m | sampling step when two anchors are compared |
| engine | SIGMA_AUTO / SIGMA_LINE / SIGMA_START | 0.30 / 0.60 / 10 m | observation and prior weights |
| engine | GREEN_SHARE / AMBER_SHARE | 0.85 / 0.45 | image share thresholds (calibrated so that 48A at 0.81 and 4.1 m out is not green) |
| engine | GREEN_MIN_NEIGHBOURS | 3 | every 3-neighbour parcel landed within 2.7 m; 2-neighbour ones were 4 to 9 m |
| engine | SUPPORT_MARGIN | 0.80 | runner-up pose must explain less than this share of the best |
| engine | DISTINCT_M / DISTINCT_DEG | 5 m / 3° | what counts as a different pose |
| engine | CHAIN_GATE_M / ADJUST_GATE_M | 40 m / 8 m | chain admission for ranking / for adjustment |
| engine | CHAIN_MIN_M / CHAIN_FIT_RMS | 12 m / 2.0 m | minimum shared run; fit quality for a candidate |
| engine | MAX_RESIDUAL_M | 5.0 m | above this the placement contradicts its neighbours: red |
| engine | MAX_OVERLAP_SHARE / OVERLAP_TOL_M | 0.02 / 3.0 m | pose may not sit on a neighbour, after 3 m erosion |
| engine | PRINTED_REACH_M | 10 m | a printed neighbour must be this close |
| engine | SIDE_COS | 0.38 (67.5°) | printed compass side tolerance |
| engine | LINE_REACH_M / LINE_STEP_M / LINE_DEG | 8 m / 2 m / 15° | point-to-line observations |
| engine | IMAGE_TOP_N, ADJUST_ROUNDS, WINDOW_MARGIN_M | 3, 3, 40 m | imagery on top-3 poses; solve rounds; satellite window margin |
| engine | POSE_TOL_DEG / POSE_TOL_M | 1° / 2 m | candidate de-duplication |
| gcp | CORNER_TURN_MIN / MATCH_DIST / TEAM_SIGMA | 20° / 1.5 m / 0.20 m | what is a corner; matched to an edge; weight of a team GCP |
| match | TURN_TOL | 10° | matched corners must bend the same way |
| match (args) | tol_abs, tol_rel, min_len, pass_deg | 0.30 m, 1 %, 8 m, 5° | chainage matching tolerances |
| raster | zoom / resolution / READ_ATTEMPTS / MAX_EMPTY_SHARE | 20 / 0.15 m / 3 / 0.02 | export window |
| edges (args) | min_len_m | 3.0 m | LSD segments shorter than this dropped |
| topology | MIN_AREA / GUARD_FRACTION / CLEAN_TOL / MIN_WIDTH / gap_tol | 0.02 m² / 0.25 / 0.02 m / 0.20 m / 1.20 m | noise floor; refuse edits eating a plot; cleaning |
| window | BUFFER_M, COARSE_STEP, FINE_STEP, SHIFT_STEP, SHIFT_RANGE, TWIN_RATIO, AREA_GATE | 40 m, 2°, 0.5°, 2 m, 6 m, 1.15, 3.0 | Puvi start-pose search (reference only now) |
| evaluate | REPORT_ONLY_RMS | 3.0 m | manual parcels looser than this are reported, not counted |
| fmb_to_geojson | A0_H, PT_MM, size_groups(min_count=4) | 3370, 25.4/72, 4 | drawing-size heuristics; glyph template scores are per entry in `glyphlib.json` (213 templates with per-template thresholds, e.g. 0.58, 0.571) |

Assumptions baked in, beyond the numbers: rigid placement only (scale exactly 1); sheets drawn north-up so printed sides rotate with the pose; the team's `.points` file is the best pose of a hand-placed parcel; a railway strip keeps its shape; Puvi is never a source of position; GCP corners must turn at least 20°; satellite edges are boundaries (they are often roads, see 47B).

### 1.3 Where it fails, as measured

- **One neighbour cannot orient a parcel.** 47B from 48A alone: right edge choice only after the printed sides were used; rotation off until the strip line constraint was added.
- **Long straight boundaries slide.** 42A: two poses explained 348 m and 347 m of boundary.
- **Image evidence prefers roads.** 47B: satellite share 0.47 at the wrong pose, 0.12 at the team's. Observable on 3 of 15 parcels only.
- **The team's references are affine.** 47B is stretched 13 to 16 %; anchor rms up to 3.05 m; 42B and 43B sit 20 m apart though the sheets name each other.
- **Transcription is manual.** 18 of 19 villages have no `nb_override`, and the tool cannot start without it.
- **Sheet conversion fails on fragments.** 22 fragment conversions and 22 Puvi identity conflicts in the corridor audit.
- Leave-one-out with 14 neighbours: 10 of 15 within 3 m. From a single seed: 171 at 3.2 m, 47B disputed.

### 1.4 What should stay rule-based, and why

Keep deterministic: **C** (anchor pose from `.points`: a closed-form fit), **D** (chain matching is exact geometry), **G** (least-squares adjustment: explainable, sigma per parcel), **I** (topology), **K** (review). Replacing these with a model would lose exactness and explainability without gaining accuracy; the sheet lengths are the ground truth and the adjustment already respects them.

Candidates for ML: **A** (glyph reading), **B** (transcription of numbers and sides), **F** (boundary evidence from imagery), **H** (GCP corner localisation), **J** (verdict). **E** could take a learned ranker later, once J's features are trusted.

## 2. Step 2: hybrid architecture, stage by stage

The six stages in your brief mapped onto this pipeline. Compute is stated for CPU on this machine unless noted; there is no GPU here.

### 2.1 Image-to-reference matching (your stage 1) = our F

What "reference" means here is important: the input is a **vector sheet**, not a scanned raster, so classic image-to-image matching (SIFT/ORB, SuperPoint+LightGlue, LoFTR) has nothing to match on the sheet side. The evidence problem is **boundary-in-imagery**, so this stage collapses into 2.4 (segmentation) plus geometric alignment.

| Candidate | Pros | Cons / failure modes | Data | Compute |
|---|---|---|---|---|
| a. LSD segments + share score (today) | zero data, fast (9 s per window), explainable | finds roads and building edges, not plot boundaries; observable on 3/15 parcels | none | CPU |
| b. Boundary probability raster from a segmentation model (2.4), then the same share/alignment on it | uses real boundary evidence; the alignment code stays | needs a model; hedges, walls and dry-season fields look alike | labels from your placements | CPU inference feasible at 0.15 m/px on 300 m windows (seconds with a small U-Net) |
| c. Learned image-to-vector matching (render the sheet, match to imagery with LoFTR-style features) | end-to-end | no pretrained model for sketch-to-satellite; would need thousands of pairs; opaque | thousands of pairs you do not have | GPU |

**Recommendation: b.** Keep `align.py`'s share, margin and observability tests, but feed them a boundary-probability raster instead of LSD segments. Training data: your placed parcels rasterised as boundary masks over the exported windows (see Step 3); pretrained backbone, fine-tuned.

### 2.2 Control-feature detection (your stage 2) = our H

FMB sheets have no grid intersections or benchmarks; the control features are **parcel corners** on the ground.

| Candidate | Pros | Cons | Data | Compute |
|---|---|---|---|---|
| a. Corner = sheet corner nearest an LSD segment (today) | free | matches roads; 1.5 m tolerance is a guess | none | CPU |
| b. Corner localisation on the boundary-probability raster: for each sheet corner, search a small window for the local maximum of the corner response (Harris on the boundary map) | reuses 2.1 b; sub-metre when the boundary is visible | fails under trees, at compound walls | none beyond 2.4 | CPU |
| c. Object detector (YOLO / RT-DETR) trained on "parcel corner" boxes | direct | corners are not objects with appearance; boxes are arbitrary; needs thousands of labelled corners | thousands of labels | GPU to train |

**Recommendation: b**, and keep writing the `.points` file exactly as today so the team's Georeferencer workflow is unchanged. No detector.

### 2.3 Text and label reading (your stage 3) = our A and B

| Candidate | Pros | Cons | Data | Compute |
|---|---|---|---|---|
| a. Glyph templates (today, 213 templates) | fast, offline | brittle: "17" for 170, village labels read as surveys; per-template thresholds | maintained by hand | CPU |
| b. PaddleOCR (PP-OCRv4, pretrained) on rendered sheet crops around each label position the PDF already gives, plus a rule layer (normalise `S.No`, reject `V.No`) | strong on printed digits, zero training to start; the PDF gives exact text positions so detection is free | Tamil script in labels; needs a confidence gate; 20 to 50 ms per crop on CPU | none to start; ~200 labelled crops to fine-tune if needed | CPU |
| c. TrOCR / a vision-language model per sheet (numbers + sides in one pass) | reads sides and neighbours together, tolerant of layout | slower, needs an API or GPU; harder to audit | the two-reader JSON for Kizhikaranai is already a labelled set | GPU or API |

**Recommendation: b for numbers and dimensions (stage A); c for the neighbour-and-side transcription (stage B)**, run twice with two readers and kept only where both agree, exactly like the manual workflow that produced `nb_override_35_04_077.json`. The Kizhikaranai JSON (15 sheets, 60+ entries) is the validation set.

### 2.4 Boundary and feature extraction (your stage 4) = new, feeds 2.1 and 2.2

| Candidate | Pros | Cons | Data | Compute |
|---|---|---|---|---|
| a. LSD (today) | none beyond zero cost | see 2.1 a | none | CPU |
| b. U-Net / DeepLabV3+ with a pretrained ImageNet encoder, fine-tuned on boundary masks from your placed parcels; output = boundary probability per pixel | small, fast on CPU, trains on hundreds of windows, directly the evidence we need | learns the corridor's look (dry fields, walls); seasonal drift; needs re-labelling per region | 15 Kizhikaranai + 11 Kolathur parcels now; every village you place adds more | train on GPU (cloud or a colleague's card), infer on CPU |
| c. SAM / SAM 2 zero-shot with the sheet polygon as a prompt (box or points) | no training; surprisingly good on field-like regions | not boundary-aware in the cadastral sense; heavy on CPU (seconds per prompt, ViT-B) | none | CPU slow, GPU fine |

**Recommendation: b, with c as a zero-shot first check before any training** (one afternoon: prompt SAM with the top-3 poses' polygons and see whether its masks separate the poses better than LSD does). If SAM already discriminates, we may defer training.

### 2.5 Transformation estimation (your stage 5) = our C, D, G

| Candidate | Pros | Cons |
|---|---|---|
| a. Rigid pose per sheet + block adjustment with Huber loss (today) | keeps every FMB length; explainable sigmas; already handles fixed and free parcels and line constraints | cannot express a genuinely distorted sheet |
| b. RANSAC/MAGSAC affine or polynomial per sheet from image-derived GCPs | standard georeferencing | throws away the dimension guarantee; with 4 to 8 corner GCPs the model is under-determined and noisy |
| c. Thin-plate spline | flexible | worst of all for cadastral: silently bends boundaries |

**Recommendation: a, unchanged.** MAGSAC is still worth adding **inside** the GCP step as an outlier filter over image-derived GCPs before they enter the adjustment (cheap, no data).

### 2.6 Quality control and confidence (your stage 6) = our J

| Candidate | Pros | Cons | Data |
|---|---|---|---|
| a. Hand rule (today) | explainable, calibrated on 15 parcels | brittle; two greens out of 15; will not transfer to other villages unchanged | 15 rows |
| b. Gradient-boosted classifier (scikit-learn, already installed) on the row features (support, runner-up ratio, residual, neighbours, anchor partners, side agreement, image share, in-window, pass number) predicting P(error < 3 m), with isotonic calibration | small data is fine; gives a probability, not a colour; feature importance is readable | needs a few hundred labelled rows; must be split by village | every leave-one-out and seed run you judge produces rows |
| c. Conformal prediction on top of b | guaranteed coverage of the auto-accept set | same data needs | same |

**Recommendation: b now, c once there are 300+ labelled rows.** Labels come from your verdicts and from the measured error against parcels you place afterwards. The colour rule stays as the fallback and as the explanation shown next to the probability.

### 2.7 Summary: one recommendation per stage and what it needs

| Stage | Recommendation | Pretrained enough? | Data you would need |
|---|---|---|---|
| Sheet numbers, dimensions | PaddleOCR + rules | yes to start | ~200 crops to fine-tune if digit confusions persist |
| Neighbour and side transcription | vision model, two readers, agreement only | yes | Kizhikaranai JSON as validation; spot-check 1 village |
| Boundary evidence | SAM zero-shot check, then fine-tuned U-Net | SAM yes; U-Net needs fine-tuning | boundary masks from your placed parcels: 26 now, hundreds after rollout |
| GCP corners | corner response on the boundary map | n/a | none |
| Placement, adjustment, topology, review | keep deterministic | n/a | none |
| Confidence | gradient boosting + calibration | n/a | 300+ judged rows, split by village |

Incremental order, one PR each: (1) confidence classifier (no new dependency, immediate value), (2) SAM zero-shot boundary check, (3) OCR for stage A, (4) transcription reader for stage B, (5) U-Net if SAM is not enough.

## 3. Clarifying questions before Step 3

1. **Accuracy target.** Your own placements disagree by 1 to 3 m. Is "within 3 m of a hand placement" the acceptance, or do you want a stricter target measured against independent check points (a GNSS survey of a few corners)? Without independent points, the ceiling is the anchors' own consistency.
2. **Reference for stretched parcels.** For 47B the FMB-exact placement disagrees with your stretched one. Which is truth for training and scoring: the sheet's dimensions, or the hand placement as saved?
3. **Compute.** There is no GPU here and no PyTorch. Training (U-Net, any fine-tuning) needs a GPU: a cloud instance for a few hours, or a machine in the office? CPU-only inference is fine for everything recommended.
4. **Imagery rights.** Google Satellite tiles are used as an internal review basemap. Training a model on them and keeping windows on disk is a licensing question for FarmwiseAI; do you have a licensed imagery source (Bhuvan, Maxar, drone orthos) we should use instead?
5. **Transcription scope.** Should stage B read only what the sheet prints (numbers, sides), or also plot numbers and dimensions inside the parcel (which stage A already reads from the PDF text)?
6. **Acceptance flow.** Do you want the confidence score to drive auto-accept (parcels above P = 0.9 written without your review) or to order the review queue only?
7. **Kolathur.** Its 11 manual parcels have no neighbour transcription. May I run the two-reader transcription there so the training and validation set is two villages, not one?

## 4. What I did not verify

- PaddleOCR, SAM and PyTorch are not installed; their CPU timings above are from their published benchmarks, not measured here.
- The corridor's other villages' sheet quality (fragments, scale mix) beyond the 22 fragment conversions already counted.
