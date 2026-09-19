"""Shared-boundary matching between two sheet outlines, and the observations it yields.

The chainage walk is ported from the 2026-09-17 prototype: walk both rings by running distance from
every start-corner pair, match corners at equal running distance, pass through a corner that turns
less than `pass_deg` (the two sheets draw a near-straight corner differently), and end the walk
where the boundaries really diverge. Every run is returned, not only the longest, because returning
only the longest hid true matches.
"""
from dataclasses import dataclass

import numpy as np
from shapely.geometry import LineString, Point

from . import sheets


@dataclass
class Chain:
    length: float
    pairs: list
    n_pairs: int
    corners: int


@dataclass
class PairObs:
    """Two sheet vertices that are the same ground point."""
    a: str
    b: str
    pa: tuple
    pb: tuple
    sigma: float


@dataclass
class LineObs:
    """A sheet boundary sample that must lie on a fixed anchor's boundary line."""
    survey: str
    p_local: tuple
    line: tuple
    sigma: float
    distance_now: float


TURN_TOL = 10.0     # degrees: matched corners must bend the same way, or the rings have parted


def common_chains(vq, vp, tol_abs=0.30, tol_rel=0.01, min_len=8.0, pass_deg=5.0, turn_tol=TURN_TOL):
    """Every run of consecutive equal-length edges shared by rings vq and vp, longest first.

    The walk ends at a matched vertex where the two rings turn differently: that is the corner
    where the shared boundary stops and each parcel goes its own way. Without this check the
    2026-09-19 run walked 46A's outline two edges past its boundary with 43B, because the next
    edges happened to be 28 m and 27 m long, and the pairs it added were 52 m and 137 m apart.
    """
    def tol(x):
        return tol_abs + tol_rel * max(x, 1.0)

    eq = sheets.edge_lengths(vq)
    tq = sheets.turn_angles(vq)
    sq = sheets.signed_turns(vq)
    nq = len(eq)
    seen, out = set(), []
    for reverse in (True, False):
        vp2 = vp[::-1] if reverse else vp
        ep = sheets.edge_lengths(vp2)
        tp = sheets.turn_angles(vp2)
        sp = sheets.signed_turns(vp2)
        np_len = len(ep)
        for i in range(nq):
            for j in range(np_len):
                pairs = [(vq[i], vp2[j])]
                cq = cp = matched = 0.0
                iq = ip = corners = 0
                while iq < nq and ip < np_len:
                    q_next = cq + eq[(i + iq) % nq]
                    p_next = cp + ep[(j + ip) % np_len]
                    if abs(q_next - p_next) <= tol(max(q_next, p_next)):
                        cq, cp, iq, ip = q_next, p_next, iq + 1, ip + 1
                        pairs.append((vq[(i + iq) % nq], vp2[(j + ip) % np_len]))
                        matched = min(cq, cp)
                        kq, kp = (i + iq) % nq, (j + ip) % np_len
                        if tq[kq] >= pass_deg or tp[kp] >= pass_deg:
                            if abs(sq[kq] - sp[kp]) > turn_tol:
                                break               # a shared corner, then the rings part
                            corners += 1
                        if iq >= nq or ip >= np_len:
                            break
                    elif q_next < p_next:
                        iq += 1
                        cq = q_next
                        if tq[(i + iq) % nq] > pass_deg:
                            break
                    else:
                        ip += 1
                        cp = p_next
                        if tp[(j + ip) % np_len] > pass_deg:
                            break
                if len(pairs) >= 2 and matched >= min_len:
                    key = tuple(sorted((round(a[0], 2), round(a[1], 2), round(b[0], 2), round(b[1], 2))
                                       for a, b in pairs))
                    if key not in seen:
                        seen.add(key)
                        out.append(Chain(matched, pairs, len(pairs), corners))
    out.sort(key=lambda c: -c.length)
    return out


def chain_observations(chain, a, b, sigma):
    return [PairObs(a, b, tuple(pa), tuple(pb), sigma) for pa, pb in chain.pairs]


def line_observations(survey, samples, anchor_boundary, sigma, max_dist=2.0):
    """Point-to-line observations against a fixed anchor's boundary.

    Point-to-line rather than point-to-point because the anchor's own geometry is affine: its
    vertices sit 1.4-6.4 m from any rigid placement of its sheet, but the line it draws is still the
    boundary. Only samples already within max_dist are used, so a wrong pose cannot recruit the far
    side of the parcel.
    """
    geoms = anchor_boundary.geoms if anchor_boundary.geom_type.startswith("Multi") else [anchor_boundary]
    out = []
    for p in np.asarray(samples, float):
        pt = Point(p)
        best = None
        for g in geoms:
            coords = list(g.coords)
            for k in range(len(coords) - 1):
                d = LineString([coords[k], coords[k + 1]]).distance(pt)
                if best is None or d < best[0]:
                    best = (d, (coords[k], coords[k + 1]))
        if best is None or best[0] > max_dist:
            continue
        out.append(LineObs(survey, (float(p[0]), float(p[1])), best[1], sigma, float(best[0])))
    return out
