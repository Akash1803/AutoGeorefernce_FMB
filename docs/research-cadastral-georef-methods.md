# Research note: georeferencing cadastral sheets, rigid pose plus block adjustment versus the alternatives

Date: 2026-09-21. Purpose: a go/no-go on keeping the rigid-pose-plus-block-adjustment core
before building the evaluation harness and a review-ordering confidence model on it.

Method: 101 agents in a deep-research workflow (5 search angles, 19 sources fetched, 94 claims
extracted, 25 verified by three independent adversarial votes each: 19 confirmed, 6 refuted).
Only confirmed claims appear below; refuted ones are listed at the end so they are not cited by
mistake. Claims without a measured figure are marked unquantified.

## Verdict

**Go, read as "no evidence to change the core", not as "settled".** The published record
supports least-squares block adjustment of adjoining cadastral sheets with shared-boundary tie
points as an established method. No alternative core (affine or polynomial fits,
rubber-sheeting, image-to-map feature matching, learned boundary extraction) has been shown to
reach cadastral accuracy in comparable conditions; that is different from being shown not to,
and the learned methods in particular were tested by nobody on line-drawing sheets like these.
Three conditions:

1. Every accuracy figure is reported as fit-to-reference with the accuracy-cap sentence.
   Already in PR 0.
2. Hold out whole villages. Already in PR 0.
3. Test the rigid assumption instead of assuming it: no source validates scale fixed at 1 as the
   final placement, so the harness should run a similarity (free-scale) arm beside the rigid
   arm and report the fitted scale distribution. As a measurement only; outputs stay scale 1.
   Proposed as a small addition to PR 0 (below), for approval.

## What the evidence says

**Block adjustment (high confidence).** Shmutter and Doytsher (1992) treat each digitised
cadastral map as a photograph and the neighbours as a photogrammetric block; points on shared
boundaries are tie points, traverse points are control. Klebanov and Doytsher (2009) implement
it with a four-parameter similarity per block: on a synthetic 300-block array with 61 control
points, inter-block MSE fell from 0.51 m to 0.08 m and error against injected truth was 0.10 to
0.11 m versus 0.29 to 0.38 m for the practised averaging method. On 11 real 1:2500 blocks the
companion FIG 2009 paper reports MSE 1.71 m for per-block-then-average against 0.90 m
(similarity) and 0.61 m (affine) for simultaneous adjustment; the reference is internal
residuals plus re-surveyed points, not a GNSS check set. The gain shrinks as control thins:
0.14 versus 0.44 m with 23 control points per 100 parcels, 0.32 versus 0.42 m with 6. The
zero-surveyed-control regime, which is ours, is untested in any paper. [1] [2]

**Rigid placement (high confidence, and an absence of evidence).** Rigid translation plus
rotation appears in the literature only as assembly or initialisation: the Como 1:2000 mosaics
(roto-translated, then a global polynomial fit) [3], the ISRO-SAC Kanchipuram village sheets
(mosaicked by translation and rotation, then 40 tie points to a DGPS-rectified IKONOS image, no
residual reported: unquantified) [4], and the Cimahi pipeline (rigid ICP first, similarity
last) [5]. No source measures scale fixed at 1 as the final transformation, so the rigid choice
is neither supported nor refuted. The sources that free scale do so for plane-table sheets with
paper deformation; whether chain-and-offset FMB sheets with printed dimensions are scale-true
is unquantified everywhere. One team concedes that mosaicking before fitting cost accuracy
compared with per-sheet georeferencing [3].

**Local warping is rejected (high confidence).** Thin-plate spline showed significant local
deformations where control was sparse [3] (qualitative, unquantified); rubber-sheeting at 2 m
and 3 m link distances produced unwanted local deformations and overlaps and did not preserve
topology [6]; a post-hoc affine with 30 control points reached 1.23 to 1.51 m at the fitting
points but made only minor changes to parcel boundaries and left most slivers [6]. Locally
adaptive adjustment helps only with dense verified control (623 points in the Slovenian case),
which we do not have.

**Achieved accuracy without GNSS check points (high confidence).** The reported figures cluster
at about 1 mm at map scale, which is 0.5 to 2 m at FMB scales: Czech 1:2880 sheets 3.1 to
4.3 m against the current cadastre [6]; Como 1:2000 sheets 5 to 10 m against municipal
cartography [3]; Cimahi 0.4 to 0.5 m against the same SAM boundaries that drove the adjustment
[5]. Every figure is bounded by its reference. The 1 to 3 m disagreement between our hand
placements sits in this band, and the accuracy of Google Satellite tiles in rural Tamil Nadu is
unquantified in every source.

**Indian sources (high confidence that they report nothing usable).** The NRDMS/KSCST village
information system guideline (2016 draft) prescribes RMS below 1 m only with dual-frequency GPS
at boundary stones and reports no achieved value [7]. The ISRO-SAC/Anna University Kanchipuram
paper (2006) gives DGPS-referenced image accuracy (0.56/0.67 m at GCPs, 1.56/1.62 m at check
points, while its conclusion states under 0.5 m) and no positional figure for the placed
cadastral mosaic; it validates lengths and areas only [4]. The Joniganuru paper (2015, a
journal on Beall's list) has no residuals and an arithmetically inconsistent table [8].
Nothing specific to Tamil Nadu FMB, Puvi, Survey of India or DILRMP accuracy surfaced, which
may reflect search limits rather than absence.

**Feature matching and learned extraction (medium confidence).** Demonstrated only as coarse
seeding or as a reference to fit to: SuperPoint plus SuperGlue on historical Jerusalem maps
reaches about 1 percent of the map diagonal, tens of metres [9]; automated matching on a
Bosnian 1:6250 sheet yields initial points but cannot find all correspondences (abstract only,
unquantified) [10]; SAM on 5 cm UAV orthophotos served as the reference, not the placement
engine [5]. Nothing shows cadastral-grade placement of line-drawing sheets against satellite
imagery, and nothing was run on a laptop CPU.

**Confidence scores ordering review (medium confidence, one source).** In [9] the number of
surviving keypoint matches predicted accuracy monotonically and a threshold accepted 76 percent
of maps while rejecting content-poor ones; thresholds were chosen in-sample, so any FMB model
must be validated on held-out villages. This supports a model that orders review rather than
replacing it, which is decision 6.

## What this changes for us

- No evidence to change the core. Rigid pose, shared-boundary observations and block
  adjustment against hand-placed anchors is the closest published analogue (Technion lineage),
  applied in the one regime nobody has measured: no surveyed control at all. The question stays
  open; the harness below is what would answer it.
- Add a scale-free comparison arm to the harness, measurement only: for every placed parcel,
  also fit a similarity transform to the same observations and log the fitted scale and the
  residual change. If the scales cluster tightly around 1, the rigid assumption is confirmed on
  our sheets; if they do not, the sheets are not scale-true and that is a finding, not a licence
  to rescale. Two new row columns (`sim_scale`, `sim_residual_m`), one function, one test.
  Awaiting approval.
- The accuracy cap is real and known: with no independent points the pipeline can prove
  internal consistency, not absolute accuracy. A handful of RTK or DGPS points at boundary
  stones, or Survey of India or DILRMP control, would lift it; that is a decision for the team,
  not the tool.
- The block-adjustment gain decays with distance from control [1]. Our pass number already
  records distance from the seed; the harness should keep reporting error against pass number,
  which the ring runs already show (3 m one hop out, 4 to 6 m two hops, 13 m three hops).

## Caveats from the verification

No confirmed source matches the FMB situation directly. All accuracy figures are map-to-map or
fit-to-reference. Klebanov's evidence is synthetic with Gaussian noise; the claim that block
adjustment gives no absolute-accuracy advantage without control was voted 1 to 2 and is not
confirmed, only the measured thinning trend is. Two sources were verified from abstracts only
[1] [10]. Some contradiction searches ran with the search budget exhausted. Learned matching is
a fast-moving area; the tens-of-metre regime may improve with fine-tuning that none of the
papers attempted, and that we cannot do on this laptop.

## Open questions the harness can answer

1. Are the portal's FMB rasters scale-true? The similarity arm above measures it.
2. What is the absolute accuracy of the hand placements and of Google tiles here? Needs a few
   independent points; unknown until then.
3. How does the adjustment behave when anchors disagree by 1 to 3 m among themselves: converge
   to their mean, inherit, or amplify at ring edges? Leave-one-village-out and the existing
   leave-one-out runs are the first measurement anywhere.

## Refuted claims, not to be cited

- That Cimahi's rigid stage contributed only 2 to 4 percent of the improvement and per-parcel
  similarity the bulk (0 to 3).
- That Cimahi's pipeline is a rigid-pose-plus-block-adjustment design close to ours (0 to 3).
- The specific 3 to 4 times improvement figure attributed to Klebanov's Table 3 (0 to 3).
- That block adjustment gives no absolute-accuracy advantage without independent control (1 to 2).
- That the ISRO-SAC authors identify mosaicking as the major difficulty (0 to 3).
- That the Bosnian result directly limits feature matching as a replacement core (0 to 3).

## Sources

1. Shmutter, Doytsher. Geomatica 46(3), 1992. https://cdnsciencepub.com/doi/10.1139/geomat-1992-0029 (abstract only)
2. Klebanov, Doytsher. Nordic Journal of Surveying and Real Estate Research 4, 2009. https://journal.fi/njs/article/view/2554
3. Brovelli, Minghini. e-Perimetron 7(3), 2012. http://www.e-perimetron.org/Vol_7_3/Brovelli_Minghini.pdf
4. Jayaprasad et al. ISPRS Archives XXXVI-4, Goa 2006. https://www.isprs.org/proceedings/xxxvi/part4/wg-iv-9-5.pdf
5. Suwardhi et al. ISPRS Archives XLVIII-2/W11, 2025. https://isprs-archives.copernicus.org/articles/XLVIII-2-W11-2025/277/2025/isprs-archives-XLVIII-2-W11-2025-277-2025.pdf
6. Kratochvilova, Cajthaml. Scientific Reports 15:26994, 2025. https://www.nature.com/articles/s41598-025-12235-9
7. NRDMS/KSCST village information system guidelines, draft 2016. https://www.kscst.org.in/vis_files/Cadastral%20_mapping_VIS_guidelines.pdf
8. Padma et al. IJESRT, July 2015 (low reliability). https://www.ijesrt.com/
9. Vaienti, di Lenardo, Kaplan. Cartography and Geographic Information Science, 2025. https://www.tandfonline.com/doi/full/10.1080/15230406.2025.2566789
10. Tuno et al. Springer LNNS, IAT 2024. https://link.springer.com/chapter/10.1007/978-3-031-71694-2_27 (abstract only)

Also consulted for the accuracy-without-GNSS angle: Bulletin of Geodetic Sciences (SciELO), ISPRS Archives XLIII-B3-2020 p.1333, PMC3791001, Survey Review 2015, Remote Sensing 14(16):4086, and the KSCST base-map SOP.
