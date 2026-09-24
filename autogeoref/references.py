"""What a placed parcel is measured against: the team's own placement, read two ways.

Reference A, FMB-exact: the sheet at the rigid pose fitted to the team's ``.points`` file, or to
the hand geometry when that fit is worse than the rms limit. Every printed length survives.
Reference B, as saved: the hand geometry exactly as the team saved it, affine and vertex-edited.

The two disagree wherever the team stretched the sheet. That disagreement is measured here in
metres along the outline, and it decides whether the parcel can carry a label (Step 3, section 3):
a label means "within the acceptance limit of the reference", which is undefined when the two
readings of the reference differ by more than that limit.
"""
import csv
import logging
import math
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

from . import anchors, config as configmod, fit, paths, sheets, visible

log = logging.getLogger(__name__)

STRETCHED_COLUMNS = ["village", "survey", "hand_file", "source", "rms_m", "area_scale",
                     "stretch_long_pct", "stretch_short_pct", "extent_m", "centroid_gap_m",
                     "disagree_mean_m", "disagree_p95_m", "disagree_max_m",
                     "disagree_corner_mean_m", "disagree_corner_max_m", "corners_matched", "band", "rule", "reason"]
CORNER_REACH_M = 10.0    # an FMB corner with no hand vertex within this counts as moved this far
DISPUTED_COLUMNS = ["village", "survey", "reason"]
REFERENCE_COLUMNS = STRETCHED_COLUMNS[:-2] + ["stretched", "stretched_reason", "disputed", "disputed_reason"]


@dataclass
class Reference:
    village: str
    survey: str
    hand_file: str
    source: str                    # "points" or "geometry"
    theta: float
    tx: float
    ty: float
    rms_m: float
    n_points: int
    area_scale: float              # sqrt(hand area / sheet area)
    stretch_long_pct: float        # minimum rotated rectangle, long side, hand vs FMB-exact
    stretch_short_pct: float
    extent_m: float                # long side of the hand geometry's rectangle
    centroid_gap_m: float
    disagree_mean_m: float         # FMB-exact outline sampled every metre, distance to the hand outline
    disagree_p95_m: float
    disagree_max_m: float
    disagree_corner_mean_m: float  # FMB-exact corners to the nearest hand vertex: sees a slide along an edge
    disagree_corner_max_m: float
    corners_matched: int
    band: str                      # "<=5%", "5-10%", ">10%" by area scale
    stretched: bool
    stretched_reason: str
    disputed: bool
    disputed_reason: str
    geom_fmb: BaseGeometry = field(repr=False, compare=False, default=None)
    geom_hand: BaseGeometry = field(repr=False, compare=False, default=None)

    def as_row(self) -> Dict[str, object]:
        return {f.name: getattr(self, f.name) for f in fields(self) if not f.name.startswith("geom_")}


def _largest(geom: BaseGeometry) -> BaseGeometry:
    if geom.geom_type == "MultiPolygon":
        return max(geom.geoms, key=lambda p: p.area)
    return geom


def _rectangle_sides(geom: BaseGeometry) -> Tuple[float, float]:
    xs = list(_largest(geom).minimum_rotated_rectangle.exterior.coords)
    sides = sorted(Point(xs[i]).distance(Point(xs[i + 1])) for i in range(4))
    return sides[-1], sides[0]


def stretch_rule(cfg: configmod.Config, area_scale: float, stretch_long_pct: float,
                 stretch_short_pct: float, rms_m: float, disagree_mean_m: float,
                 disagree_corner_mean_m: float = 0.0) -> Tuple[bool, str]:
    """Apply the configured rule. Returns (stretched, reason).

    The displacement rule takes the larger of two means: outline samples to the hand outline
    (blind to a slide along a straight edge) and FMB corners to the nearest hand vertex (sees it).
    """
    st = cfg.stretch
    if rms_m > st.rms_m:
        return True, "rigid fit rms %.2f m above %.1f m" % (rms_m, st.rms_m)
    if st.rule == "pct":
        off = abs(area_scale - 1.0) * 100
        if off > st.pct:
            return True, "area scale %.3f is %.1f %% off (limit %.0f %%)" % (area_scale, off, st.pct)
    elif st.rule == "axis":
        off = max(abs(stretch_long_pct), abs(stretch_short_pct))
        if off > st.pct:
            return True, "rectangle side off by %.1f %% (limit %.0f %%)" % (off, st.pct)
    elif st.rule == "displacement":
        worst = max(disagree_mean_m, disagree_corner_mean_m)
        if worst > st.displacement_mean_m:
            return True, "FMB-exact and hand placement disagree by %.2f m mean (outline %.2f, corners %.2f; limit %.1f m)" % (
                worst, disagree_mean_m, disagree_corner_mean_m, st.displacement_mean_m)
    else:
        raise ValueError("unknown stretch rule %r" % st.rule)
    return False, ""


def band_of(area_scale: float, band_pct: float, pct: float) -> str:
    off = abs(area_scale - 1.0) * 100
    if off <= band_pct:
        return "<=%.0f%%" % band_pct
    if off <= pct:
        return "%.0f-%.0f%%" % (band_pct, pct)
    return ">%.0f%%" % pct


def disputed_reasons(village: str, survey: str) -> List[str]:
    """Why a survey cannot be a reference at all, from the corridor audit and the run reports."""
    reasons: List[str] = []
    qc = paths.PROJECT / "FMB_Vector" / "puvi_vs_sheet_qc.csv"
    if qc.exists():
        with qc.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("village_code") == village and r.get("survey_no") == survey \
                        and "identity mismatch" in (r.get("verdict") or ""):
                    reasons.append("Puvi and the portal describe different parcels")
    rc = paths.vector_dir(village) / "rail_conflicts.csv"
    if rc.exists():
        with rc.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if survey in (r.get("survey_a"), r.get("survey_b")):
                    reasons.append("railway land conflict with %s" % (r.get("survey_b") if r.get("survey_a") == survey else r.get("survey_a")))
    manual = paths.vector_dir(village) / "disputed.csv"
    if manual.exists():
        with manual.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("survey") == survey:
                    reasons.append(r.get("reason") or "listed in disputed.csv")
    corrections = paths.vector_dir(village) / "sheet_corrections.csv"
    if corrections.exists():
        with corrections.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("survey") == survey and (r.get("edge") or "") != "plots":
                    reasons.append("outline corrected by hand (sheet_corrections.csv, verified_by=%s)"
                                   % (r.get("verified_by") or "unverified"))
    return reasons


def reference_for(village: str, survey: str, cfg: Optional[configmod.Config] = None) -> Optional[Reference]:
    """Both references for one hand-placed survey, or None when the team has not placed it."""
    cfg = cfg or configmod.load()
    hand = paths.manual_files(village, survey)
    if not hand:
        return None
    pose = anchors.pose_from_points(village, survey)
    source = "points"
    if pose is None or pose[2] > cfg.stretch.rms_m:
        pose = anchors.pose_from_geometry(village, survey, hand[0])
        source = "geometry"
    theta, t, rms, _scale, _aniso, n = pose
    sheet = sheets.load_sheet(village, survey)
    geom_fmb = sheets.dissolve([(None, fit.apply_pose(g, theta, t)) for _p, g in sheet])
    hand_gdf = visible.read_hand(hand[0])
    if hand_gdf.crs is not None and hand_gdf.crs.to_epsg() != 32644:
        hand_gdf = hand_gdf.to_crs(32644)
    geom_hand = unary_union([g.buffer(0) for g in hand_gdf.geometry if g is not None and not g.is_empty])
    sheet_area = sum(g.area for _p, g in sheet)
    area_scale = math.sqrt(geom_hand.area / sheet_area) if sheet_area > 0 else float("nan")
    long_a, short_a = _rectangle_sides(geom_fmb)
    long_b, short_b = _rectangle_sides(geom_hand)
    stretch_long = (long_b / long_a - 1.0) * 100 if long_a > 0 else float("nan")
    stretch_short = (short_b / short_a - 1.0) * 100 if short_a > 0 else float("nan")
    samples, _bearings = sheets.sample_outline(sheets.outline(sheet), 1.0)
    ground = fit.transform_points(samples, theta, t)
    boundary = geom_hand.boundary
    d = np.array([boundary.distance(Point(p)) for p in ground]) if len(ground) else np.zeros(1)
    hand_vertices = [Point(c) for poly in (geom_hand.geoms if geom_hand.geom_type == "MultiPolygon" else [geom_hand])
                     for c in poly.exterior.coords[:-1]]
    corner_d = []
    for c in geom_fmb.exterior.coords[:-1]:
        pc = Point(c)
        corner_d.append(min([pc.distance(q) for q in hand_vertices] + [CORNER_REACH_M]))
    corner_d = np.array(corner_d) if corner_d else np.zeros(1)
    stretched, reason = stretch_rule(cfg, area_scale, stretch_long, stretch_short, float(rms), float(d.mean()),
                                     float(corner_d.mean()))
    reasons = disputed_reasons(village, survey)
    ref = Reference(
        village=village, survey=survey, hand_file=str(hand[0]), source=source,
        theta=float(theta), tx=float(t[0]), ty=float(t[1]), rms_m=round(float(rms), 3), n_points=int(n),
        area_scale=round(area_scale, 4), stretch_long_pct=round(stretch_long, 2),
        stretch_short_pct=round(stretch_short, 2), extent_m=round(long_b, 2),
        centroid_gap_m=round(float(geom_fmb.centroid.distance(geom_hand.centroid)), 3),
        disagree_mean_m=round(float(d.mean()), 3), disagree_p95_m=round(float(np.percentile(d, 95)), 3),
        disagree_max_m=round(float(d.max()), 3),
        disagree_corner_mean_m=round(float(corner_d.mean()), 3), disagree_corner_max_m=round(float(corner_d.max()), 3),
        corners_matched=int((corner_d < CORNER_REACH_M).sum()),
        band=band_of(area_scale, cfg.stretch.band_pct, cfg.stretch.pct),
        stretched=bool(stretched), stretched_reason=reason,
        disputed=bool(reasons), disputed_reason="; ".join(reasons),
        geom_fmb=geom_fmb, geom_hand=geom_hand)
    log.debug("reference %s/%s: %s rms %.2f scale %.3f disagree %.2f m stretched=%s",
              village, survey, source, rms, area_scale, d.mean(), stretched)
    return ref


def references(village: str, cfg: Optional[configmod.Config] = None) -> Dict[str, Reference]:
    """Every hand-placed survey of the village, keyed by survey."""
    cfg = cfg or configmod.load()
    out: Dict[str, Reference] = {}
    for survey in paths.surveys_with_sheets(village):
        if not paths.manual_files(village, survey):
            continue
        try:
            ref = reference_for(village, survey, cfg)
        except Exception as exc:                       # a broken hand file must not stop the run
            log.warning("reference for %s/%s failed: %s", village, survey, exc)
            continue
        if ref is not None:
            out[survey] = ref
    return out


def _merge_csv(path: Path, columns: List[str], village: str, new_rows: List[Dict[str, object]]) -> Path:
    """Replace this village's rows in a corridor-wide CSV, keeping the other villages' rows."""
    kept: List[Dict[str, object]] = []
    if path.exists():
        with path.open(encoding="utf-8") as fh:
            kept = [r for r in csv.DictReader(fh) if r.get("village") != village]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for r in kept + new_rows:
            w.writerow({c: r.get(c, "") for c in columns})
    return path


def write_exclusions(village: str, refs: Dict[str, Reference], cfg: Optional[configmod.Config] = None,
                     out_dir: Optional[Path] = None) -> Tuple[Path, Path, Path]:
    """excluded_stretched.csv, excluded_disputed.csv and references.csv under _logs/eval.

    All three are corridor-wide files; this village's rows are replaced, other villages' kept.
    """
    cfg = cfg or configmod.load()
    out_dir = out_dir or (paths.logs_dir() / "eval")
    stretched_rows = [dict(r.as_row(), rule=cfg.stretch.rule, reason=r.stretched_reason)
                      for r in refs.values() if r.stretched]
    disputed_rows = [{"village": village, "survey": r.survey, "reason": r.disputed_reason}
                     for r in refs.values() if r.disputed]
    all_rows = [r.as_row() for r in refs.values()]
    p1 = _merge_csv(out_dir / "excluded_stretched.csv", STRETCHED_COLUMNS, village, stretched_rows)
    p2 = _merge_csv(out_dir / "excluded_disputed.csv", DISPUTED_COLUMNS, village, disputed_rows)
    p3 = _merge_csv(out_dir / "references.csv", REFERENCE_COLUMNS, village, all_rows)
    log.info("%s: %d references, %d stretched, %d disputed", village, len(refs), len(stretched_rows), len(disputed_rows))
    return p1, p2, p3
