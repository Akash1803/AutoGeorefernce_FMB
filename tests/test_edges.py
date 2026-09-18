import math

import numpy as np
import pytest

from autogeoref import edges

SPIKE = r"D:\code\FMB_to_GeoJSON\autogeoref\scratch_2026-09-17\spike_gcp\kizhi_sat_z20.tif"


def test_detect_finds_thousands_of_segments_on_the_spike_raster():
    segs = edges.detect(SPIKE, min_len_m=4.0)
    assert 1500 < len(segs) < 6000, "the 2026-09-18 spike found 2643 segments >= 4 m"
    assert all(s.length >= 4.0 for s in segs)
    assert all(0 <= s.bearing < math.pi for s in segs)


def test_segments_carry_ground_coordinates():
    segs = edges.detect(SPIKE, min_len_m=4.0)
    xs = [s.p0[0] for s in segs]
    ys = [s.p0[1] for s in segs]
    assert 392000 < min(xs) and max(xs) < 394000
    assert 1412000 < min(ys) and max(ys) < 1415000


def test_index_returns_distance_and_bearing_of_the_nearest_segment():
    idx = edges.SegmentIndex([edges.Segment((0.0, 0.0), (100.0, 0.0), 0.0, 100.0)], step=0.5)
    d, b = idx.query(np.array([[50.0, 2.0], [50.0, 40.0]]))
    assert d[0] == pytest.approx(2.0, abs=0.3)
    assert b[0] == pytest.approx(0.0, abs=1e-6)
    assert not np.isfinite(d[1]), "beyond max_dist there is no nearest segment"


def test_directions_groups_parallel_segments():
    segs = [edges.Segment((0, 0), (100, 0), 0.0, 100.0),
            edges.Segment((0, 10), (60, 10), 0.0, 60.0),
            edges.Segment((0, 0), (0, 40), math.pi / 2, 40.0)]
    got = edges.directions(segs, tol_deg=10.0)
    assert len(got) == 2
    assert got[0][1] == pytest.approx(160.0)


def test_bearings_are_maths_angles_not_compass_bearings():
    """s.bearing is measured anticlockwise from +x. The corridor runs at compass 25.8 deg."""
    east = edges.Segment((0, 0), (10, 0), 0.0, 10.0)
    assert edges.compass(east.bearing) == pytest.approx(90.0)
    north = edges.Segment((0, 0), (0, 10), math.pi / 2, 10.0)
    assert edges.compass(north.bearing) == pytest.approx(0.0)


def test_the_corridor_parallel_bias_is_real():
    """Imagery here is dominated by the corridor: it fixes rotation and cross-track, not along-track."""
    segs = edges.detect(SPIKE, min_len_m=4.0)
    total = sum(s.length for s in segs)
    corridor = math.radians((90 - 25.8) % 180)          # compass 25.8 deg as a maths angle
    par = sum(s.length for s in segs
              if min(abs(s.bearing - corridor), math.pi - abs(s.bearing - corridor)) <= math.radians(15))
    assert par / total > 0.6, "the 2026-09-18 spike measured 72 %"
    groups = edges.directions(segs, tol_deg=10.0)
    assert edges.compass(groups[0][0]) == pytest.approx(25.8, abs=3.0), "the corridor dominates"
    assert groups[0][1] / total > 0.5


def test_detect_can_be_limited_to_a_window():
    whole = edges.detect(SPIKE, min_len_m=4.0)
    part = edges.detect(SPIKE, bounds=(392800, 1413200, 392950, 1413350), min_len_m=4.0)
    assert 0 < len(part) < len(whole)
    for s in part:
        assert 392700 < s.p0[0] < 393050 and 1413100 < s.p0[1] < 1413450
