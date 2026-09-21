"""One feature row per parcel per run, and the village-level split.

Rows are what PR 3 may learn from and what every report is computed from. They hold the
placement features the engine already computes, the rule-based verdict, the reference the parcel
is measured against (see ``references``), the errors against both readings of that reference, and
the label. Nothing derived from imagery is stored except scalars and the window's hash.

The split is by village and only by village: features of one parcel are computed from its
neighbours, and every run of a village re-uses the same geometry, so rows of one village are never
divided between train and test.
"""
import csv
import datetime
import logging
import re
import subprocess
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from shapely.geometry.base import BaseGeometry

from . import config as configmod, paths, references as refmod, review, sheets

log = logging.getLogger(__name__)

ROWS_FILE = "rows.csv"


@dataclass
class FeatureRow:
    # identity
    run_id: str
    run_at: str
    tool_git_sha: str
    village: str
    survey: str
    mode: str                      # loo, seed, full
    seeds: str
    pass_no: Optional[int]
    imagery_source: str
    window_sha256: str
    backfilled: int
    # placement features
    method: str
    n_candidates: Optional[int]
    support_m: Optional[float]
    support_next_m: Optional[float]
    runner_up_ratio: Optional[float]
    boundary_rms_m: Optional[float]
    n_neighbours: Optional[int]
    anchor_partners: Optional[int]
    side_ok: Optional[int]
    side_bad: Optional[int]
    contradictions: Optional[int]
    printed_far_n: Optional[int]
    share: Optional[float]
    margin: Optional[float]
    observable: Optional[int]
    ambiguous: Optional[int]
    in_window: Optional[int]
    sigma_pos_m: Optional[float]
    sigma_head_deg: Optional[float]
    heading_deg: Optional[float]
    pose_tx: Optional[float]
    pose_ty: Optional[float]
    topo_max_move_m: Optional[float]
    topo_clipped: Optional[int]
    # rule-based verdict
    colour: str
    confidence_rule: Optional[int]
    # reference
    ref_source: str
    ref_rms_m: Optional[float]
    ref_area_scale: Optional[float]
    ref_stretch_long_pct: Optional[float]
    ref_stretch_short_pct: Optional[float]
    ref_disagree_mean_m: Optional[float]
    ref_disagree_p95_m: Optional[float]
    ref_disagree_corner_mean_m: Optional[float]
    stretch_band: str
    stretched: Optional[int]
    stretched_reason: str
    disputed: Optional[int]
    disputed_reason: str
    # errors
    err_fmb_m: Optional[float]
    err_hand_m: Optional[float]
    heading_err_deg: Optional[float]
    max_vertex_err_m: Optional[float]
    # label
    within_3m: Optional[int]
    notes: str

    def as_dict(self) -> Dict[str, object]:
        return {f.name: ("" if getattr(self, f.name) is None else getattr(self, f.name)) for f in fields(self)}


FIELDS: List[str] = [f.name for f in fields(FeatureRow)]


def _num(x: object) -> Optional[float]:
    if x is None or x == "":
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _int(x: object) -> Optional[int]:
    v = _num(x)
    return None if v is None else int(round(v))


def _flag(x: object) -> Optional[int]:
    if x is None or x == "":
        return None
    if isinstance(x, str):
        return 1 if x.strip().lower() in ("true", "1", "yes") else 0
    return 1 if x else 0


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=str(paths.CODE),
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return ""


def _topo(notes: str) -> Tuple[Optional[float], Optional[int]]:
    m = re.search(r"conform max ([0-9.]+) m", notes or "")
    return (float(m.group(1)) if m else None), (1 if "clipped" in (notes or "") else 0)


def heading_diff(a: float, b: float) -> float:
    return abs(((a - b) + 180.0) % 360.0 - 180.0)


def errors_against(placed: Optional[BaseGeometry], heading: Optional[float], ref: Optional[refmod.Reference]
                   ) -> Dict[str, Optional[float]]:
    """Centroid, heading and largest-vertex errors of a placed body against both references."""
    out: Dict[str, Optional[float]] = {"err_fmb_m": None, "err_hand_m": None, "heading_err_deg": None,
                                       "max_vertex_err_m": None}
    if placed is None or placed.is_empty or ref is None:
        return out
    out["err_fmb_m"] = round(float(placed.centroid.distance(ref.geom_fmb.centroid)), 3)
    out["err_hand_m"] = round(float(placed.centroid.distance(ref.geom_hand.centroid)), 3)
    out["max_vertex_err_m"] = round(float(placed.hausdorff_distance(ref.geom_fmb)), 3)
    if heading is not None:
        out["heading_err_deg"] = round(heading_diff(heading, ref.theta), 3)
    return out


def label_for(err_fmb_m: Optional[float], ref: Optional[refmod.Reference], cfg: configmod.Config) -> Optional[int]:
    """1 or 0 when the parcel has a usable reference, None when it cannot carry a label."""
    if ref is None or ref.stretched or ref.disputed or err_fmb_m is None:
        return None
    return 1 if err_fmb_m <= cfg.acceptance_m else 0


def row_from_status(village: str, run_id: str, run_at: str, mode: str, seeds: Iterable[str],
                    status: Dict[str, str], ref: Optional[refmod.Reference], cfg: configmod.Config,
                    placed: Optional[BaseGeometry], imagery_source: str = "", window_sha256: str = "",
                    tool_git_sha: str = "", backfilled: int = 0) -> FeatureRow:
    """One row from one status-CSV line plus the parcel's written geometry."""
    notes = status.get("notes") or ""
    support = _num(status.get("support_m"))
    nxt = _num(status.get("support_next_m"))
    heading = _num(status.get("heading_deg"))
    errs = errors_against(placed, heading, ref)
    topo_move, topo_clip = _topo(notes)
    return FeatureRow(
        run_id=run_id, run_at=run_at, tool_git_sha=tool_git_sha, village=village,
        survey=str(status.get("survey")), mode=mode, seeds="+".join(str(s) for s in seeds),
        pass_no=_int(status.get("pass")), imagery_source=imagery_source, window_sha256=window_sha256,
        backfilled=int(backfilled),
        method=status.get("method") or "", n_candidates=_int(status.get("n_candidates")),
        support_m=support, support_next_m=nxt,
        runner_up_ratio=(round(nxt / support, 4) if support and nxt is not None else None),
        boundary_rms_m=_num(status.get("boundary_rms_m")),
        n_neighbours=(len([x for x in (status.get("neighbours") or "").split(",") if x])
                      if status.get("neighbours") is not None else None),
        anchor_partners=_int(status.get("anchor_partners")), side_ok=_int(status.get("side_ok")),
        side_bad=_int(status.get("side_bad")), contradictions=_int(status.get("contradictions")),
        printed_far_n=len([x for x in (status.get("printed_far") or "").split(",") if x]),
        share=_num(status.get("share")), margin=_num(status.get("margin")),
        observable=_flag(status.get("observable")), ambiguous=_flag(status.get("ambiguous")),
        in_window=_flag(status.get("in_window")), sigma_pos_m=_num(status.get("sigma_pos_m")),
        sigma_head_deg=_num(status.get("sigma_head_deg")), heading_deg=heading,
        pose_tx=_num(status.get("pose_tx")), pose_ty=_num(status.get("pose_ty")),
        topo_max_move_m=topo_move, topo_clipped=topo_clip,
        colour=status.get("colour") or "", confidence_rule=_int(status.get("confidence")),
        ref_source=(ref.source if ref else "none"), ref_rms_m=(ref.rms_m if ref else None),
        ref_area_scale=(ref.area_scale if ref else None),
        ref_stretch_long_pct=(ref.stretch_long_pct if ref else None),
        ref_stretch_short_pct=(ref.stretch_short_pct if ref else None),
        ref_disagree_mean_m=(ref.disagree_mean_m if ref else None),
        ref_disagree_p95_m=(ref.disagree_p95_m if ref else None),
        ref_disagree_corner_mean_m=(ref.disagree_corner_mean_m if ref else None),
        stretch_band=(ref.band if ref else ""), stretched=(int(ref.stretched) if ref else None),
        stretched_reason=(ref.stretched_reason if ref else ""),
        disputed=(int(ref.disputed) if ref else None), disputed_reason=(ref.disputed_reason if ref else ""),
        err_fmb_m=errs["err_fmb_m"], err_hand_m=errs["err_hand_m"],
        heading_err_deg=errs["heading_err_deg"], max_vertex_err_m=errs["max_vertex_err_m"],
        within_3m=label_for(errs["err_fmb_m"], ref, cfg), notes=notes)


def placed_body(village: str, survey: str) -> Optional[BaseGeometry]:
    """The dissolved body of the tool's written file for a survey, under the current project root."""
    out = paths.output_path(village, survey)
    if not out.exists():
        return None
    import geopandas as gpd
    got = gpd.read_file(out, layer="parcels")
    return sheets.dissolve([(None, g) for g in got.geometry])


def rows_from_run(village: str, run_id: str, mode: str, seeds: Iterable[str], cfg: configmod.Config,
                  refs: Dict[str, refmod.Reference], run_at: Optional[str] = None,
                  imagery: Optional[Dict[str, str]] = None, tool_git_sha: Optional[str] = None) -> List[FeatureRow]:
    """Rows for every placed parcel of the run whose status CSV is under the current project root.

    ``refs`` come from the real project (the team's files); the caller may have redirected
    ``paths.PROJECT`` to a work copy for the status and outputs, which is the point.
    """
    run_at = run_at or datetime.datetime.now().isoformat(timespec="seconds")
    sha = git_sha() if tool_git_sha is None else tool_git_sha
    imagery = imagery or {}
    rows: List[FeatureRow] = []
    for status in review.read_status(village):
        if status.get("status") != "placed":
            continue
        survey = status["survey"]
        ref = refs.get(survey)
        placed = placed_body(village, survey)
        rows.append(row_from_status(village, run_id, run_at, mode, seeds, status, ref, cfg, placed,
                                    imagery_source=imagery.get("source", ""),
                                    window_sha256=imagery.get("sha256", ""), tool_git_sha=sha))
    log.info("%s run %s: %d rows, %d labelled", village, run_id, len(rows),
             sum(1 for r in rows if r.within_3m is not None))
    return rows


def rows_path(path: Optional[Path] = None) -> Path:
    return Path(path) if path else paths.logs_dir() / "eval" / ROWS_FILE


def append_rows(rows: Iterable[FeatureRow], path: Optional[Path] = None) -> Path:
    p = rows_path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    new = not p.exists() or p.stat().st_size == 0
    with p.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        for r in rows:
            w.writerow(r.as_dict())
    return p


def load_rows(path: Optional[Path] = None) -> List[Dict[str, str]]:
    p = rows_path(path)
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def is_labelled(row: Dict[str, object]) -> bool:
    return str(row.get("within_3m", "")) in ("0", "1")


def split_by_village(rows: List[Dict[str, object]], test_village: str
                     ) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """(train, test): every row of ``test_village`` is test, every other row is train."""
    test = [r for r in rows if r.get("village") == test_village]
    train = [r for r in rows if r.get("village") != test_village]
    return train, test


def folds(rows: List[Dict[str, object]]) -> List[Tuple[str, List[Dict[str, object]], List[Dict[str, object]]]]:
    """Leave-one-village-out over the villages that have at least one labelled row."""
    villages = sorted({str(r.get("village")) for r in rows if is_labelled(r)})
    out = []
    for v in villages:
        train, test = split_by_village(rows, v)
        out.append((v, train, test))
    return out


def backfill_leave_one_out(csv_path: Path, village: str, run_id: str, cfg: configmod.Config,
                           refs: Dict[str, refmod.Reference]) -> List[FeatureRow]:
    """Rows from a leave-one-out CSV written before this harness existed; missing fields stay empty."""
    rows: List[FeatureRow] = []
    with Path(csv_path).open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            survey = r.get("survey")
            ref = refs.get(survey)
            status = {"survey": survey, "colour": r.get("colour"), "method": r.get("method"),
                      "share": r.get("share"), "observable": r.get("observable"),
                      "boundary_rms_m": r.get("boundary_rms_m"), "neighbours": r.get("neighbours"),
                      "margin": r.get("margin"), "notes": "backfilled from %s" % Path(csv_path).name}
            row = row_from_status(village, run_id, "", "loo", [], status, ref, cfg, None, backfilled=1)
            err = _num(r.get("centroid_error_m"))
            row.err_fmb_m = err
            row.heading_err_deg = _num(r.get("heading_error_deg"))
            row.within_3m = label_for(err, ref, cfg)
            rows.append(row)
    return rows
