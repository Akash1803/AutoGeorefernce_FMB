import geopandas as gpd
import numpy as np
from shapely import affinity
from shapely.geometry import Polygon

from autogeoref import fmb_stage2


def sq(x0, y0, w=30.0, h=30.0):
    return Polygon([(x0, y0), (x0 + w, y0), (x0 + w, y0 + h), (x0, y0 + h)])


def test_neighbour_pairs_come_from_the_base_fabric():
    base = gpd.GeoDataFrame({"survey_no": ["1", "2", "3"]},
                            geometry=[sq(0, 0), sq(30.5, 0), sq(200, 0)], crs=32644)
    assert fmb_stage2.neighbour_pairs(base) == [("1", "2")]


def test_two_free_sheets_close_their_shared_edge():
    # B is placed 2.5 m east of where it meets A; the shared edge should close
    geoms = {"A": sq(0, 0), "B": sq(32.5, 0)}
    deltas = fmb_stage2.adjust(geoms, [("A", "B")], free={"A", "B"})
    a = fmb_stage2.apply_delta(geoms["A"], deltas["A"])
    b = fmb_stage2.apply_delta(geoms["B"], deltas["B"])
    assert a.exterior.distance(b.exterior) < 1.0
    assert all(d["flag"] == "" for d in deltas.values())


def test_a_fixed_neighbour_never_moves():
    geoms = {"A": sq(0, 0), "B": sq(32.5, 0)}
    deltas = fmb_stage2.adjust(geoms, [("A", "B")], free={"B"})
    assert "A" not in deltas, "only free parcels get a delta"
    b = fmb_stage2.apply_delta(geoms["B"], deltas["B"])
    # B walked west towards the fixed A
    assert min(x for x, _y in b.exterior.coords) < 31.5


def test_a_parcel_with_no_nearby_neighbour_stays_put():
    geoms = {"A": sq(0, 0), "B": sq(80, 0)}      # 50 m apart: a real disagreement, not a seam
    deltas = fmb_stage2.adjust(geoms, [("A", "B")], free={"A", "B"})
    assert deltas["B"]["moved_m"] < 0.01 and deltas["A"]["moved_m"] < 0.01


def test_the_adjustment_is_rigid():
    geoms = {"A": sq(0, 0), "B": affinity.rotate(sq(31.5, 0), 3.0, origin="centroid")}
    deltas = fmb_stage2.adjust(geoms, [("A", "B")], free={"B"})
    b = fmb_stage2.apply_delta(geoms["B"], deltas["B"])
    assert abs(b.area - geoms["B"].area) < 1e-6
    assert abs(b.length - geoms["B"].length) < 1e-6


def test_a_move_beyond_the_guard_is_refused():
    # closing this seam would take a 7.9 m move; that is not a seam, and the guard says so
    geoms = {"A": sq(0, 0), "B": sq(37.9, 0)}
    deltas = fmb_stage2.adjust(geoms, [("A", "B")], free={"B"})
    assert deltas["B"]["flag"].startswith("refused")
    assert deltas["B"]["moved_m"] == 0.0
    b = fmb_stage2.apply_delta(geoms["B"], deltas["B"])
    assert b.equals_exact(geoms["B"], 1e-9), "a refused move leaves the parcel where stage 1 put it"


def test_the_gate_only_adopts_an_improvement():
    better = [{"before_m": 3.0, "after_m": 2.0}, {"before_m": 4.0, "after_m": 3.5},
              {"before_m": 2.0, "after_m": 2.0}]
    worse = [{"before_m": 3.0, "after_m": 5.0}, {"before_m": 4.0, "after_m": 4.1},
             {"before_m": 2.0, "after_m": 2.0}]
    assert fmb_stage2.gate_verdict(better)["adopt"] is True
    assert fmb_stage2.gate_verdict(worse)["adopt"] is False
    assert fmb_stage2.gate_verdict([])["adopt"] is False, "no evidence is not a pass"


def test_the_thresholds_are_the_ones_the_sop_states():
    assert fmb_stage2.NEAR_M == 8.0
    assert fmb_stage2.PRIOR_POS_M == 3.0
    assert fmb_stage2.MAX_MOVE_M == 6.0
