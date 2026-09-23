"""Shape QC of converted sheets, run before any placement.

One drawn segment 9.9 m short of its printed length (613) moved four parcels on 2026-09-21, and
offset arrows polygonised as plot lines cut slivers off real plots (47A). This module flags what
can be seen in the converted geometry alone: slivers, spikes, unnumbered plots, invalid or
multi-part outlines, plots that overlap, holes between plots, and thin lines that cross the
survey boundary away from their endpoints (the arrow signature). Printed-length checks need the
dimension labels and are a separate step.
"""
import csv
import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional

from shapely.geometry import Point, shape
from shapely.ops import unary_union

from . import paths

log = logging.getLogger(__name__)

SLIVER_M2 = 3.0
THIN_COMPACTNESS = 0.08
SPIKE_DEG = 8.0
SPIKE_ARM_M = 3.0
QC_COLUMNS = ["village", "survey", "plots", "outline_m2", "n_issues", "issues"]


def _compactness(g) -> float:
    return 4 * math.pi * g.area / g.length ** 2 if g.length else 0.0


def check_plots(plots: List[tuple]) -> List[str]:
    """Issues in a list of (plot_no, polygon) pairs of one sheet."""
    issues: List[str] = []
    for name, g in plots:
        if not g.is_valid:
            issues.append("plot %s invalid" % name)
        if g.area < SLIVER_M2:
            issues.append("plot %s sliver %.1f m2" % (name, g.area))
        elif _compactness(g) < THIN_COMPACTNESS:
            issues.append("plot %s very thin (compactness %.2f, %.0f m2)" % (name, _compactness(g), g.area))
        if name in ("None", "nan", "", None):
            issues.append("unnumbered plot %.1f m2" % g.area)
        ring = list(g.exterior.coords)
        n = len(ring) - 1
        for i in range(n):
            a, b, c = ring[i - 1], ring[i], ring[(i + 1) % n]
            v1 = (a[0] - b[0], a[1] - b[1])
            v2 = (c[0] - b[0], c[1] - b[1])
            l1, l2 = math.hypot(*v1), math.hypot(*v2)
            if min(l1, l2) < SPIKE_ARM_M:
                continue
            ang = math.degrees(math.acos(max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1]) / (l1 * l2)))))
            if ang < SPIKE_DEG:
                issues.append("plot %s spike %.0f deg, arms %.1f/%.1f m" % (name, ang, l1, l2))
    for i in range(len(plots)):
        for j in range(i + 1, len(plots)):
            ov = plots[i][1].intersection(plots[j][1]).area
            if ov > 0.5:
                issues.append("plots %s and %s overlap %.1f m2" % (plots[i][0], plots[j][0], ov))
    if plots:
        body = unary_union([g for _, g in plots])
        parts = list(body.geoms) if body.geom_type == "MultiPolygon" else [body]
        holes = sum(len(p.interiors) for p in parts)
        if holes:
            issues.append("%d hole(s) between plots" % holes)
        if len(parts) > 1:
            issues.append("outline in %d parts" % len(parts))
    return issues


def check_lines(lines: List[dict]) -> List[str]:
    """Thin lines crossing the survey boundary away from their endpoints: the arrow signature."""
    issues: List[str] = []
    bnd = [shape(f["geometry"]) for f in lines if f["properties"].get("layer") == "survey_boundary_main"]
    for f in lines:
        if f["properties"].get("layer") != "subdivision_line":
            continue
        g = shape(f["geometry"])
        for b in bnd:
            x = g.intersection(b)
            pts = [x] if x.geom_type == "Point" else list(getattr(x, "geoms", []))
            for pt in pts:
                if pt.geom_type != "Point":
                    continue
                ends = (Point(g.coords[0]), Point(g.coords[-1]), Point(b.coords[0]), Point(b.coords[-1]))
                if min(pt.distance(e) for e in ends) > 0.05:
                    issues.append("line %.1f m crosses the boundary at (%.1f, %.1f)" % (g.length, pt.x, pt.y))
    return issues


def check_sheet(village: str, survey: str) -> Dict[str, object]:
    d = json.loads((paths.vector_dir(village) / ("%s_parcels.geojson" % survey)).read_text(encoding="utf-8"))
    plots = [(str(f["properties"].get("plot_no")), shape(f["geometry"])) for f in d["features"]]
    issues = check_plots(plots)
    lp = paths.vector_dir(village) / ("%s_lines.geojson" % survey)
    if lp.exists():
        issues += check_lines(json.loads(lp.read_text(encoding="utf-8"))["features"])
    outline = unary_union([g for _, g in plots]).area if plots else 0.0
    return {"village": village, "survey": survey, "plots": len(plots), "outline_m2": round(outline, 1),
            "n_issues": len(issues), "issues": issues}


def qc_village(village: str, surveys: Optional[List[str]] = None, out_csv: Optional[Path] = None) -> Dict[str, object]:
    surveys = list(surveys) if surveys else paths.surveys_with_sheets(village)
    rows = [check_sheet(village, s) for s in surveys]
    p = Path(out_csv) if out_csv else paths.vector_dir(village) / "sheet_qc.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=QC_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, issues="; ".join(r["issues"])))
    log.info("%s: %d sheets, %d with issues -> %s", village, len(rows), sum(1 for r in rows if r["issues"]), p)
    return {"village": village, "sheets": rows, "csv": str(p)}
