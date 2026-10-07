"""Place the next ring of FMB sheets onto Akash's placed parcels (the seeds), the way the approved
Thirukatchur batches of 25-28 Sep 2026 were placed:

* seeds are his hand-placed / approved layers in QGIS group Manual_Georeferenced/<Village>_Village,
  exactly as QGIS shows them; they are never changed;
* a sheet joins the ring when it and a placed parcel print each other's number (printed_numbers);
* the sheet keeps its PDF shape (rotation + shift only); its edges facing a placed neighbour are laid
  on that neighbour's boundary - candidate poses from runs of equal edges (match.common_chains) and
  from edge-on-edge alignment, anchored at a corner;
* a pose is judged by where the sheet's own printed numbers land: a number that lands inside a placed
  parcel with another number is a contradiction, so is a long shared boundary with a parcel neither
  sheet names; shared boundary with named neighbours is the support;
* a sheet is placed only when one pose clearly wins; otherwise it is asked, with the reason.
"""
import math
import re
from dataclasses import dataclass, field

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import shapely.ops
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import unary_union

from . import fit, match, paths, printed_numbers, sheets

UTM = 32644
CHAIN_TOL_ABS, CHAIN_TOL_REL, CHAIN_MIN_M = 0.5, 0.02, 5.0
EDGE_MIN_M = 5.0            # shortest edge used for edge-on-edge poses
SUPPORT_M = 2.0             # boundary within this of a neighbour's boundary is shared boundary (hand
                            # placements differ from their own sheets by 1-4 m: 68B rms 2.5 m)
SAMPLE_M = 0.5
ERODE_M = 2.0               # overlap is measured between parcels shrunk by this (hand placements disagree by 1-4 m)
MAX_OVERLAP = 0.03          # share of the sheet's area
CONTRA_M = 8.0              # shared boundary with a parcel neither sheet names, beyond this, contradicts
LABEL_NEAR_M = 25.0         # a named neighbour's printed number must land within this (+ its distance
                            # from the sheet's outline) of that neighbour
LABEL_DEPTH_M = 10.0        # a number of an unplaced survey may land this deep inside a placed parcel
MIN_SUPPORT_M = 15.0
RUNNER_UP = 0.80            # the best pose must beat every distinct pose by this ratio
CLUSTER_M = 4.0             # poses whose corners all land within this are one answer (hand work differs 1-4 m)
CROSS_DEG = 25.0            # support in two directions this far apart pins the position along an edge
REFINE_ROUNDS = 3
STRIP_ELONGATION = 5.0      # length / width beyond which a sheet is treated as a strip (rail, road)
ROT_WINDOW_DEG = 25.0        # a pose may turn this much (+ the seeds' own spread) from the seeds' heading
REFINE_TOP = 8               # distinct poses refined before the best is chosen
JUNCTION_REACH_M = 6.0      # a sheet corner this close to a junction of two placed parcels is pinned to it
REFINE_REACH_M = 2.5        # outline samples this close to a named neighbour pull onto its line


# --------------------------------------------------------------------------- data
@dataclass
class Sheet:
    survey: str
    plots: list                 # [(properties, Polygon in sheet metres)]
    outline: list               # clockwise ring, sheet metres
    polygon: Polygon
    labels: list                # [printed_numbers.Label]
    meta: dict = field(default_factory=dict)


@dataclass
class Decision:
    survey: str
    status: str                 # placed | ask
    reason: str = ""
    theta: float = 0.0
    t: tuple = (0.0, 0.0)
    metrics: dict = field(default_factory=dict)
    options: list = field(default_factory=list)     # [(theta, t, shared boundary m, how found)], best first


def clockwise(coords):
    pts = [tuple(map(float, p[:2])) for p in coords]
    if len(pts) > 1 and math.dist(pts[0], pts[-1]) < 1e-9:
        pts = pts[:-1]
    out = [pts[0]]
    for p in pts[1:]:
        if math.dist(p, out[-1]) > 1e-6:
            out.append(p)
    return out[::-1] if Polygon(out).exterior.is_ccw else out


def largest(g):
    if g is None or g.is_empty:
        return None
    if isinstance(g, Polygon):
        return g
    parts = [p for p in getattr(g, "geoms", []) if isinstance(p, Polygon)]
    return max(parts, key=lambda p: p.area) if parts else None


def load_sheet(village, survey, table):
    plots = [(props, poly) for props, poly in sheets.load_sheet(village, survey) if poly.area > 0]
    import json
    meta = json.loads(paths.sheet_path(village, survey).read_text(encoding="utf-8")).get("metadata", {})
    ring = sheets.outline(plots)
    return Sheet(str(survey), plots, ring, Polygon(ring), table.get(str(survey), []), meta)


def pose_points(P, theta, t):
    return fit.transform_points(P, theta, t)


def ground_polygon(sheet, theta, t):
    return Polygon(pose_points(sheet.outline, theta, t))


# --------------------------------------------------------------------------- context
class Context:
    """Placed parcels (seeds + this batch so far) and the printed-number tables."""

    def __init__(self, placed, table, seed_ids):
        self.placed = {str(k): v for k, v in placed.items()}       # survey -> Polygon/MultiPolygon (UTM)
        self.table = table
        self.printed = printed_numbers.printed_sets(table)
        self.agreed = printed_numbers.agreed_sets(table)
        self.seed_ids = set(map(str, seed_ids))
        self.poses = {}                     # survey -> (theta, t) of its sheet on the ground, where known
        self.outlines = {}                  # survey -> sheet outline Polygon (sheet metres)
        self._rebuild()

    def _rebuild(self):
        self.ids = list(self.placed)
        self.geoms = [self.placed[s] for s in self.ids]
        self.tree = shapely.STRtree(self.geoms)
        self.union = unary_union(self.geoms)
        self.free_boundary = self.union.boundary
        self.junctions = self._junctions()

    def _junctions(self):
        """Corners on the outer boundary of the placed parcels where two of them meet: where a new sheet's
        corner must land (582's frontage corners in the Thirukatchur batches)."""
        out = []
        for i, g in enumerate(self.geoms):
            for x, y in shapely.get_coordinates(g.boundary):
                pt = Point(x, y)
                if self.free_boundary.distance(pt) > 0.2:
                    continue
                owners = {self.ids[k] for k in self.tree.query(pt.buffer(0.5))
                          if self.geoms[k].boundary.distance(pt) < 0.5}
                if len(owners) >= 2:
                    out.append((x, y))
        return np.unique(np.round(np.array(out), 2), axis=0) if out else np.zeros((0, 2))

    def add(self, survey, geom, pose=None):
        self.placed[str(survey)] = geom
        if pose is not None:
            self.poses[str(survey)] = pose
        self._rebuild()

    def rail(self):
        if not hasattr(self, "_rail"):
            try:
                r = gpd.read_file(paths.RAIL_GPKG, layer=paths.RAIL_LAYER).to_crs(UTM)
                self._rail = unary_union(r.geometry.values)
            except Exception:
                self._rail = None
        return self._rail

    def set_rotation(self, theta, spread):
        """The village's sheets are drawn near one heading (Paranur: -12..0 deg on all 12 hand-placed
        sheets); a pose turned far from it is wrong (62A at 58.9 deg, 2026-10-07)."""
        self.rotation, self.rotation_spread = theta, spread

    def rotation_ok(self, theta):
        if getattr(self, "rotation", None) is None:
            return True
        return abs(((theta - self.rotation) + 180) % 360 - 180) <= ROT_WINDOW_DEG + self.rotation_spread

    own_sheet = False       # also trust a neighbour the sheet being placed prints (70A prints 72; 72 prints nothing back)

    def named(self, a, b):
        """'agreed' | 'one-sided' | 'neither' | 'unknown' (a sheet is not read)."""
        a, b = str(a).upper(), str(b).upper()
        if a not in self.printed or b not in self.printed:
            return "unknown"
        x, y = b in self.printed[a], a in self.printed[b]
        if self.own_sheet and x and not y:
            return "agreed"
        return "agreed" if x and y else "one-sided" if x or y else "neither"


# --------------------------------------------------------------------------- poses
def candidate_poses(sheet, ctx, neighbour):
    """Rigid poses that lay part of the sheet's outline on the neighbour's free boundary."""
    N = largest(ctx.placed[neighbour])
    vq = sheet.outline
    vp = clockwise(N.exterior.coords)
    out = []
    for ch in match.common_chains(vq, vp, tol_abs=CHAIN_TOL_ABS, tol_rel=CHAIN_TOL_REL, min_len=CHAIN_MIN_M):
        P = [a for a, _b in ch.pairs]
        Q = [b for _a, b in ch.pairs]
        on_free = np.mean([ctx.free_boundary.distance(Point(q)) < 0.5 for q in Q])
        if on_free < 0.5:
            continue
        th, t, _rms, mx = fit.rigid_fit(P, Q)
        if mx <= 3.0:
            out.append((th, np.asarray(t), "chain %.0f m on %s" % (ch.length, neighbour)))
    labels = [L for L in sheet.labels if L.number == neighbour.upper()]
    reach = min((sheet.polygon.exterior.distance(Point(L.x, L.y)) for L in labels), default=None)
    for i in range(len(vq)):
        a, b = np.array(vq[i]), np.array(vq[(i + 1) % len(vq)])
        ls = float(np.linalg.norm(b - a))
        if ls < EDGE_MIN_M:
            continue
        if reach is not None:
            d = min(LineString([tuple(a), tuple(b)]).distance(Point(L.x, L.y)) for L in labels)
            if d > reach + 30.0:
                continue                       # this edge does not face the neighbour's printed number
        for j in range(len(vp)):
            p, q = np.array(vp[j]), np.array(vp[(j + 1) % len(vp)])
            lg = float(np.linalg.norm(q - p))
            if lg < EDGE_MIN_M:
                continue
            mid = Point(*(p + q) / 2)
            if ctx.free_boundary.distance(mid) > 0.5:
                continue                       # an inner seam between two placed parcels
            # two clockwise rings run along a shared edge in opposite directions: a->b along q->p
            th = math.degrees(math.atan2(*(p - q)[::-1]) - math.atan2(*(b - a)[::-1]))
            R = fit.transform_points([a, b], th, (0.0, 0.0))
            out.append((th, q - R[0], "edge %.0f m on %s (corner)" % (min(ls, lg), neighbour)))
            out.append((th, p - R[1], "edge %.0f m on %s (corner)" % (min(ls, lg), neighbour)))
    return out


def _samples(ring_pts, step=SAMPLE_M):
    pts, brg = sheets.sample_outline(ring_pts, step)
    return pts, brg


def evaluate(sheet, theta, t, ctx):
    """Support, contradictions and overlap of one pose."""
    G = ground_polygon(sheet, theta, t)
    if not G.is_valid:
        G = G.buffer(0)
    m = {"overlap": 0.0, "support": {}, "contra": {}, "label_conflicts": [], "label_far": [], "dirs": []}
    # where the sheet's printed numbers land
    for L in sheet.labels:
        if L.number == sheet.survey.upper():
            continue
        pt = Point(*pose_points([(L.x, L.y)], theta, t)[0])
        for k in ctx.tree.query(pt, predicate="within"):
            owner = ctx.ids[k]
            if owner.upper() == L.number:
                continue
            # a number printed beyond a narrow neighbour (a 17 m rail strip) lands in the parcel behind it
            named = next((x for x in ctx.ids if x.upper() == L.number), None)
            reach = sheet.polygon.exterior.distance(Point(L.x, L.y))
            if named is not None:
                if ctx.placed[named].distance(pt) <= LABEL_NEAR_M + reach:
                    continue
            elif ctx.geoms[k].boundary.distance(pt) <= LABEL_DEPTH_M + reach:
                continue        # an unplaced neighbour's number printed far out (72's "71/A" lands in 70A)
            m["label_conflicts"].append("%s lands in %s" % (L.text, owner))
    for n in ctx.ids:
        if ctx.named(sheet.survey, n) != "agreed":
            continue
        labs = [L for L in sheet.labels if L.number == n.upper()]
        if labs:
            d, reach = min((ctx.placed[n].distance(Point(*pose_points([(L.x, L.y)], theta, t)[0])),
                            sheet.polygon.exterior.distance(Point(L.x, L.y))) for L in labs)
            if d > LABEL_NEAR_M + reach:
                m["label_far"].append("%s %.0f m away" % (n, d))
    me = sheet.survey.upper()
    for n, (pth, pt_) in ctx.poses.items():
        labs = [L for L in ctx.table.get(n, []) if L.number == me]
        if not labs or n not in ctx.placed:
            continue
        reach = getattr(ctx, "_reach", {}).get(n)
        if reach is None:
            ctx._reach = getattr(ctx, "_reach", {})
            reach = ctx._reach[n] = {}
        best = None
        for L in labs:
            if (L.x, L.y) not in reach:
                outline = ctx.outlines.get(n)
                reach[(L.x, L.y)] = outline.exterior.distance(Point(L.x, L.y)) if outline is not None else 0.0
            d = G.distance(Point(*pose_points([(L.x, L.y)], pth, pt_)[0]))
            best = d - reach[(L.x, L.y)] if best is None else min(best, d - reach[(L.x, L.y)])
        if best is not None and best > LABEL_NEAR_M:
            m["label_far"].append("%s's sheet puts %s %.0f m away" % (n, sheet.survey, best + 0))
    near = [k for k in ctx.tree.query(G.buffer(SUPPORT_M + 1))]
    if near:
        U = unary_union([ctx.geoms[k] for k in near])
        inner = G.buffer(-ERODE_M)
        m["overlap"] = (inner.intersection(U.buffer(-ERODE_M)).area / G.area) if not inner.is_empty else 0.0
        pts, brg = _samples(list(G.exterior.coords)[:-1])
        P = shapely.points(pts)
        for k in near:
            n = ctx.ids[k]
            d = shapely.distance(P, ctx.geoms[k].boundary)
            hit = d < SUPPORT_M
            L = float(hit.sum()) * SAMPLE_M
            if L <= 0:
                continue
            rel = ctx.named(sheet.survey, n)
            if rel == "agreed":
                m["support"][n] = L
                m["dirs"] += list(np.degrees(brg[hit]) % 180.0)
            elif rel == "neither" and L > CONTRA_M:
                m["contra"][n] = L
    m["total"] = sum(m["support"].values())
    m["two_dirs"] = _two_directions(m["dirs"])
    m["valid"] = (m["overlap"] <= MAX_OVERLAP and not m["label_conflicts"] and not m["contra"]
                  and not m["label_far"] and m["total"] > 0)
    return m


def _two_directions(dirs, min_m=5.0):
    if not dirs:
        return False
    d = np.asarray(dirs)
    hist = {}
    for b in d:
        k = int(b // 5)
        hist[k] = hist.get(k, 0) + SAMPLE_M
    main = max(hist, key=hist.get) * 5 + 2.5
    off = np.minimum(np.abs(d - main), 180 - np.abs(d - main))
    return float((off >= CROSS_DEG).sum()) * SAMPLE_M >= min_m


def refine(sheet, theta, t, ctx, neighbours_, heading=None):
    """Rigid least squares: every outline sample near a named neighbour onto that neighbour's line."""
    for _ in range(REFINE_ROUNDS):
        G_pts, _b = _samples(pose_points(sheet.outline, theta, t), 1.0)
        local = fit.unapply_pose(G_pts, theta, t)
        obs = []
        for n in neighbours_:
            b = ctx.placed[n].boundary
            for pl, pg in zip(local, G_pts):
                pt = Point(pg)
                if b.distance(pt) > REFINE_REACH_M:
                    continue
                q = shapely.ops.nearest_points(b, pt)[0]
                seg = _segment_at(b, q)
                if seg is None:
                    continue
                obs.append(match.LineObs(sheet.survey, (float(pl[0]), float(pl[1])), seg, 0.5, float(b.distance(pt))))
        gcps = []
        if len(ctx.junctions):
            V = pose_points(sheet.outline, theta, t)
            for pl, pg in zip(sheet.outline, V):
                d = np.linalg.norm(ctx.junctions - pg, axis=1)
                k = int(np.argmin(d))
                if d[k] <= JUNCTION_REACH_M:
                    gcps.append(fit.GcpObs(sheet.survey, tuple(pl), tuple(ctx.junctions[k]), 0.7))
        if len(obs) + 2 * len(gcps) < 3:
            break
        prior = (fit.PosePrior(sheet.survey, heading, tuple(t), 50.0, 0.3) if heading is not None
                 else fit.PosePrior(sheet.survey, theta, tuple(t), 3.0, 3.0))
        sol = fit.block_adjust({sheet.survey: (theta, tuple(t))}, {}, [], obs, gcps, [prior])[sheet.survey]
        theta, t = sol["theta"], np.asarray(sol["t"])
    return theta, t


def _segment_at(boundary, q):
    """The boundary segment that point q (on the boundary) lies on."""
    lines = boundary.geoms if hasattr(boundary, "geoms") else [boundary]
    best = None
    for ln in lines:
        c = list(ln.coords)
        for k in range(len(c) - 1):
            d = LineString([c[k], c[k + 1]]).distance(q)
            if best is None or d < best[0]:
                best = (d, (c[k], c[k + 1]))
    return best[1] if best else None


def _axis_bearing(pts):
    """Direction (degrees, 0..180) of the main axis of a point set."""
    P = np.asarray(pts, float) - np.mean(pts, axis=0)
    _u, _s, vt = np.linalg.svd(P, full_matrices=False)
    return math.degrees(math.atan2(vt[0][1], vt[0][0])) % 180.0


def elongation(sheet):
    r = sheet.polygon.minimum_rotated_rectangle
    c = list(r.exterior.coords)
    a, b = math.dist(c[0], c[1]), math.dist(c[1], c[2])
    return max(a, b) / max(min(a, b), 1e-6)


def track_aligned(sheet, theta, t, ctx):
    """For a rail strip: the heading that makes the strip's long axis run along the track near it,
    or None when the track does not run through the strip (Akash: rail strips sit on the track)."""
    rail = ctx.rail()
    if rail is None or elongation(sheet) < STRIP_ELONGATION:
        return None
    G = ground_polygon(sheet, theta, t)
    inside = rail.intersection(G.buffer(10))
    if inside.is_empty or inside.length < 0.3 * math.sqrt(G.area * elongation(sheet)):
        return None
    tp = shapely.get_coordinates(inside.segmentize(5.0))
    if len(tp) < 3:
        return None
    d = (_axis_bearing(tp) - _axis_bearing(list(G.exterior.coords)[:-1]) + 90) % 180 - 90
    return theta + d


def _cluster(cands):
    """Distinct poses (best of each cluster), best first. Two poses are the same answer when no corner of
    the sheet lands more than CLUSTER_M apart: three poses of 74 (50 ha) 1-3 m apart were reported as
    'two positions' on 2026-10-07, while a 3 deg turn of a 1 km parcel is a different answer."""
    cands = sorted(cands, key=lambda c: -c["m"]["total"])
    clusters = []
    for c in cands:
        pc = pose_points(c["outline"], c["theta"], c["t"])
        for k in clusters:
            if np.max(np.linalg.norm(pc - pose_points(k["outline"], k["theta"], k["t"]), axis=1)) <= CLUSTER_M:
                break
        else:
            clusters.append(c)
    return clusters


def decide(sheet, ctx):
    """Place the sheet or say why it is asked."""
    named = [n for n in ctx.ids if ctx.named(sheet.survey, n) == "agreed"]
    if not named:
        return Decision(sheet.survey, "ask", "no placed parcel names it and is named by it")
    outline = np.asarray(sheet.outline)
    cands, seen, why = [], set(), {}
    for n in named:
        for th, t, src in candidate_poses(sheet, ctx, n):
            key = (round(th % 360, 1), round(float(t[0]), 1), round(float(t[1]), 1))
            if key in seen:
                continue
            seen.add(key)
            if not ctx.rotation_ok(th):
                continue
            m = evaluate(sheet, th, t, ctx)
            if m["valid"]:
                cands.append({"theta": th, "t": t, "m": m, "src": src, "outline": outline})
            else:
                if m["overlap"] > MAX_OVERLAP:
                    why["overlaps a placed parcel"] = why.get("overlaps a placed parcel", 0) + 1
                for x in m["label_conflicts"][:1] + m["label_far"][:1]:
                    why[x.split(" lands")[0] + (" lands in the wrong parcel" if " lands" in x else "")] = why.get(x, 0) + 1
                for x in list(m["contra"])[:1]:
                    why["runs along %s (not a neighbour)" % x] = why.get("runs along %s (not a neighbour)" % x, 0) + 1
    if not cands:
        top = sorted(why.items(), key=lambda kv: -kv[1])[:3]
        return Decision(sheet.survey, "ask", "no position fits (%s)" % "; ".join("%s x%d" % kv for kv in top),
                        metrics={"named_placed": ",".join(named)})
    # every distinct pose is pulled onto the named neighbours' lines before they are compared: a pose
    # 3 m off the right one shares little boundary until it is refined (68B, 2026-10-07)
    refined = []
    for c in _cluster(cands)[:REFINE_TOP]:
        th, t = refine(sheet, c["theta"], c["t"], ctx, named)
        m = evaluate(sheet, th, t, ctx)
        if m["valid"] and m["total"] >= c["m"]["total"] and ctx.rotation_ok(th):
            refined.append(dict(c, theta=th, t=np.asarray(t), m=m))
        else:
            refined.append(c)
    # a rail strip is turned to run along the track, then pulled onto its neighbours again
    if elongation(sheet) >= STRIP_ELONGATION:
        for k, c in enumerate(refined):
            th_t = track_aligned(sheet, c["theta"], c["t"], ctx)
            if th_t is None:
                continue
            pin = np.mean(pose_points(sheet.outline, c["theta"], c["t"]), axis=0)
            if c["m"]["support"]:                     # turn about the stretch it shares with its neighbours
                G = ground_polygon(sheet, c["theta"], c["t"])
                shared = G.exterior.intersection(unary_union([ctx.placed[n] for n in c["m"]["support"]]).buffer(SUPPORT_M))
                if not shared.is_empty:
                    pin = np.asarray(shared.centroid.coords[0])
            th2 = th_t                                   # rotate the pose about `pin`
            local_pin = fit.unapply_pose([pin], c["theta"], c["t"])[0]
            t2 = pin - fit.transform_points([local_pin], th2, (0.0, 0.0))[0]
            th3, t3 = refine(sheet, th2, t2, ctx, named, heading=th2)
            m3 = evaluate(sheet, th3, t3, ctx)
            if m3["valid"]:
                refined[k] = dict(c, theta=th3, t=np.asarray(t3), m=m3, src=c["src"] + " + along the track")
    clusters = _cluster(refined)
    best = clusters[0]
    th, t, m = best["theta"], best["t"], best["m"]
    runner = clusters[1]["m"]["total"] if len(clusters) > 1 else 0.0
    ratio = runner / m["total"] if m["total"] else 1.0
    options = [(c["theta"], tuple(c["t"]), c["m"]["total"], c["src"]) for c in clusters[:3]]
    metrics = {"support_m": round(m["total"], 1),
               "support_by": ", ".join("%s %.0f m" % (k, v) for k, v in sorted(m["support"].items(), key=lambda x: -x[1])),
               "runner_up_ratio": round(ratio, 2), "two_directions": m["two_dirs"],
               "overlap_share": round(m["overlap"], 4), "labels_far": "; ".join(m["label_far"]),
               "candidates": len(cands), "pose_from": best["src"], "named_placed": ",".join(named)}
    notes = []
    if m["total"] < MIN_SUPPORT_M:
        notes.append("only %.0f m shared with your parcels" % m["total"])
    if ratio >= RUNNER_UP:
        notes.append("a position %.0f m away fits nearly as well (%.0f vs %.0f m shared)" % (
            _pose_gap(sheet, best, clusters[1]), m["total"], runner))
    if not m["two_dirs"] and "chain" not in best["src"]:
        notes.append("touches your parcels along one straight line: its position along it is the weakest part")
    metrics["confidence"] = "check" if notes else "good"
    return Decision(sheet.survey, "placed", "; ".join(notes), th, tuple(t), metrics, options)


def _pose_gap(sheet, a, b):
    return float(np.max(np.linalg.norm(pose_points(sheet.outline, a["theta"], a["t"])
                                       - pose_points(sheet.outline, b["theta"], b["t"]), axis=1)))


# --------------------------------------------------------------------------- the ring
def next_ring(ctx, sheets_available, size):
    placed = set(ctx.placed)
    rows = []
    for s in sheets_available:
        if s in placed:
            continue
        n = [x for x in ctx.agreed.get(s, set()) if x in placed]
        if n:
            rows.append((s, len(n)))
    rows.sort(key=lambda r: (-r[1], paths.survey_sort_key(r[0])))
    return [s for s, _n in rows[:size]]


def seed_rotations(village, seeds, table):
    """{survey: heading} of each seed's sheet on the ground (rigid fit of the sheet to his outline)."""
    out = {}
    for s, g in seeds.items():
        if not paths.sheet_path(village, s).exists():
            continue
        sh = load_sheet(village, s, table)
        poly = largest(g)
        chains = match.common_chains(sh.outline, clockwise(poly.exterior.coords), tol_abs=3.0, tol_rel=0.1, min_len=5.0)
        if not chains:
            continue
        ch = chains[0]
        th, _t, rms, _mx = fit.rigid_fit([a for a, _b in ch.pairs], [b for _a, b in ch.pairs])
        if rms <= 5.0 and ch.n_pairs >= 3:
            out[s] = (((th + 180) % 360) - 180, np.asarray(_t))
    return out


def place_ring(village, seeds, table, size=8, ring=None, log=print, in_order=False, own_sheet=False):
    """Greedy: place the sheet with the clearest win, add it to the placed set, repeat. With in_order the
    sheets are placed one by one in the order given (rail strips first, Akash's rule of 2026-09-23)."""
    ctx = Context(seeds, table, seeds.keys())
    ctx.own_sheet = own_sheet
    seed_poses = seed_rotations(village, seeds, table)
    ctx.poses.update(seed_poses)
    for k in seed_poses:
        ctx.outlines[k] = load_sheet(village, k, table).polygon
    rots = {k: v[0] for k, v in seed_poses.items()}
    if len(rots) >= 2:
        v = np.array(list(rots.values()))
        med = float(np.degrees(np.arctan2(np.sin(np.radians(v)).mean(), np.cos(np.radians(v)).mean())))
        spread = float(np.max(np.abs(((v - med) + 180) % 360 - 180)))
        ctx.set_rotation(med, spread)
        log("sheet heading on the ground: %.1f deg (seeds %s), allowed %.1f +- %.1f" % (
            med, ", ".join("%s %.1f" % kv for kv in sorted(rots.items())), med, ROT_WINDOW_DEG + spread))
    available = [s for s in table if paths.sheet_path(village, s).exists()]
    ring = ring or next_ring(ctx, available, size)
    log("ring: %s" % ", ".join(ring))
    sheet_objs = {s: load_sheet(village, s, table) for s in ring}
    for k, sh_ in sheet_objs.items():
        ctx.outlines[k] = sh_.polygon
    pending, done = list(ring), {}
    if in_order:
        for s_ in ring:
            d = decide(sheet_objs[s_], ctx)
            if d.status == "placed":
                log("placed %-5s support %5.0f m (%s) runner-up %.2f from %s %s" % (
                    d.survey, d.metrics["support_m"], d.metrics["support_by"], d.metrics["runner_up_ratio"],
                    d.metrics["pose_from"], ("| " + d.reason) if d.reason else ""))
                ctx.add(d.survey, ground_polygon(sheet_objs[d.survey], d.theta, d.t), (d.theta, np.asarray(d.t)))
            done[s_] = d
        pending = []
    while pending:
        decisions = {s: decide(sheet_objs[s], ctx) for s in pending}
        ok = [d for d in decisions.values() if d.status == "placed"]
        if not ok:
            for s in pending:
                done[s] = decisions[s]
            break
        d = max(ok, key=lambda d: (d.metrics["support_m"] * (1 - d.metrics["runner_up_ratio"])))
        log("placed %-5s support %5.0f m (%s) runner-up %.2f from %s" % (
            d.survey, d.metrics["support_m"], d.metrics["support_by"], d.metrics["runner_up_ratio"], d.metrics["pose_from"]))
        ctx.add(d.survey, ground_polygon(sheet_objs[d.survey], d.theta, d.t), (d.theta, np.asarray(d.t)))
        done[d.survey] = d
        pending.remove(d.survey)
    for s in ring:
        if done[s].status == "ask":
            log("ask    %-5s %s" % (s, done[s].reason))
    return [done[s] for s in ring], sheet_objs, ctx


def plots_on_ground(sheet, theta, t):
    rows = []
    for props, poly in sheet.plots:
        g = fit.apply_pose(poly, theta, t)
        rows.append(dict(props, geometry=g))
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=UTM)
