"""Topology between parcels the tool may move. Anchors are never touched.

No snapping. The 2026-09-18 16:20 run used a 0.30 m snap across every polygon of the village; on
526A, whose 177 subdivision plots include slivers under a metre across, that collapsed two plots to
nothing and cut 153 m2 off a third. Errors are fixed only by clipping an overlap out of one side and
by giving a gap to the polygon that borders it most, and every edit is refused if it would empty a
polygon or take more than GUARD_FRACTION of its area.
"""
import csv
import itertools

from shapely.geometry import Polygon
from shapely.ops import unary_union

from . import paths

MIN_AREA = 0.02          # square metres: smaller than this is numerical noise
GUARD_FRACTION = 0.25    # refuse any edit that takes more than this share of a polygon


def _largest(geom):
    """Repair a geometry without silently discarding parts of a multipolygon."""
    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == "GeometryCollection":
        polys = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
        return unary_union(polys) if polys else Polygon()
    return geom


def check(geoms):
    keys = sorted(geoms, key=paths.survey_sort_key)
    overlaps = []
    for a, b in itertools.combinations(keys, 2):
        area = geoms[a].intersection(geoms[b]).area
        if area > MIN_AREA:
            overlaps.append([a, b, round(area, 2)])
    union = unary_union([g.buffer(0) for g in geoms.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    gaps = [round(Polygon(r).area, 2) for p in parts for r in p.interiors if Polygon(r).area > 0.01]
    return {"overlap_pairs": len(overlaps), "overlap_sqm": round(sum(o[2] for o in overlaps), 2),
            "overlaps": overlaps, "gaps": len(gaps), "gap_sqm": round(sum(gaps), 2),
            "gap_list": sorted(gaps, reverse=True)[:15],
            "invalid": [k for k, g in geoms.items() if not g.is_valid],
            "union_parts": len(parts)}


def fix(geoms, movable, rail, gap_tol=1.20):
    """Clip overlaps and fill gaps between `movable` parcels. Returns (geoms, report)."""
    out = {k: g.buffer(0) for k, g in geoms.items()}
    report = {"clips": [], "fills": [], "refused": [], "rail_conflicts": []}
    keys = sorted(out, key=paths.survey_sort_key)

    def edit(key, new_geom, why, taken):
        old = out[key]
        if new_geom is None or new_geom.is_empty:
            report["refused"].append([key, why, round(taken, 2), "would empty the polygon"])
            return False
        if old.area > 0 and taken / old.area > GUARD_FRACTION:
            report["refused"].append([key, why, round(taken, 2),
                                      "would take %.0f %% of a %.2f m2 polygon"
                                      % (100 * taken / old.area, old.area)])
            return False
        out[key] = _largest(new_geom)
        return True

    # 1. overlaps: one side yields
    for _round in range(4):
        changed = 0
        for a, b in itertools.combinations(keys, 2):
            inter = out[a].intersection(out[b])
            if inter.is_empty or inter.area <= MIN_AREA:
                continue
            ra, rb = a in rail, b in rail
            if ra != rb:
                # railway land was carved out of the older parcel, so the strip keeps its shape.
                # A parcel the tool placed is clipped to it (guarded); a team parcel is only reported.
                pair = (a, b, round(inter.area, 2))
                if pair not in report["rail_conflicts"]:
                    report["rail_conflicts"].append(pair)
                other = b if ra else a
                strip = a if ra else b
                if other in movable and edit(other, out[other].difference(out[strip]),
                                             "clip against railway land %s" % strip, inter.area):
                    report["clips"].append([other, "clipped against railway land", strip, round(inter.area, 2)])
                    changed += 1
                continue
            big, small = (a, b) if out[a].area >= out[b].area else (b, a)
            if big not in movable:
                big, small = small, big
            if big not in movable:
                continue
            if edit(big, out[big].difference(out[small]), "clip against %s" % small, inter.area):
                report["clips"].append([big, "clipped against", small, round(inter.area, 2)])
                changed += 1
        if not changed:
            break

    # 2. enclosed gaps, then thin open gaps between two surveys
    union = unary_union([g.buffer(0) for g in out.values()])
    parts = union.geoms if union.geom_type == "MultiPolygon" else [union]
    candidates = [(Polygon(r), "enclosed gap") for p in parts for r in p.interiors
                  if Polygon(r).area > MIN_AREA]
    for a, b in itertools.combinations(keys, 2):
        d = out[a].distance(out[b])
        if d <= 0 or d > gap_tol:
            continue
        strip = out[a].buffer(gap_tol).intersection(out[b].buffer(gap_tol)).difference(out[a]).difference(out[b])
        for piece in (strip.geoms if strip.geom_type.startswith("Multi") else [strip]):
            if piece.is_empty or piece.area <= MIN_AREA or piece.length == 0:
                continue
            if 2 * piece.area / piece.length > gap_tol:       # mean width: must be a sliver
                continue
            candidates.append((piece, "gap between %s and %s" % (a, b)))
    for gap, kind in candidates:
        best, best_share = None, 0.0
        for key in movable:
            share = out[key].buffer(0.05).intersection(gap).area
            if share > best_share:
                best, best_share = key, share
        if best is None:
            continue
        out[best] = _largest(unary_union([out[best], gap]))
        report["fills"].append([round(gap.area, 2), "into", best, kind])
    return out, report


def apply_to_parts(parts, fixed_body, original_body):
    """Carry a survey's topology-fixed body down to its sub-plots.

    `parts` is [(props, geom)] as placed. Every part is clipped to the fixed body (so a clipped
    overlap disappears from the plot that held it) and each filled gap piece goes to the part
    that borders it most. Returns [(props, geom, note)] with note "" / "clipped" / "filled".
    """
    out = []
    added = fixed_body.difference(original_body)
    pieces = [p for p in (added.geoms if added.geom_type.startswith("Multi") else [added])
              if not p.is_empty and p.area > MIN_AREA]
    owner = {}
    for i, piece in enumerate(pieces):
        best, best_len = None, 0.0
        for k, (_props, g) in enumerate(parts):
            shared = g.buffer(0.05).intersection(piece).area
            if shared > best_len:
                best, best_len = k, shared
        if best is not None:
            owner.setdefault(best, []).append(piece)
    for k, (props, g) in enumerate(parts):
        ng = g.intersection(fixed_body)
        note = ""
        if abs(ng.area - g.area) > MIN_AREA:
            note = "clipped"
        if k in owner:
            ng = unary_union([ng] + owner[k])
            note = (note + ",filled").strip(",")
        ng = _largest(ng)
        if ng is None or ng.is_empty or ng.geom_type not in ("Polygon", "MultiPolygon"):
            ng, note = g, "kept (edit would empty the plot)"
        out.append((props, ng, note))
    return out


def report_rail_conflicts(village, report):
    p = paths.vector_dir(village) / "rail_conflicts.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["survey_a", "survey_b", "overlap_sqm", "note"])
        for a, b, area in report.get("rail_conflicts", []):
            w.writerow([a, b, area,
                        "railway land carved out of the older parcel sheet; decide in review"])
    return p
