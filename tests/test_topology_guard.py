from shapely.geometry import Polygon

from autogeoref import topology


def _sq(x0, y0, x1, y1):
    return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def test_sanity_reports_a_parcel_that_lost_too_much_area():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20))]
    carved = [({"plot_no": "1"}, _sq(0, 0, 20, 12), "clipped")]
    parts, notes = topology.sanity(rigid, carved)
    assert notes == ["topology cut 40 %"] and parts[0][1].area == 240


def test_sanity_replaces_a_spiky_plot_with_the_rigid_plot_clipped_to_the_body():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20)), ({"plot_no": "2"}, _sq(20, 0, 30, 20))]
    # the conform step moved a needle of plot 1 into plot 2: plot 1 has a notch, plot 2 a hairpin,
    # the body (their union) is still the clean 30 x 20 rectangle
    needle = Polygon([(20, 10), (5, 10.1), (20, 10.2)])
    notched = _sq(0, 0, 20, 20).difference(needle)
    hairpin = _sq(20, 0, 30, 20).union(needle)
    parts, notes = topology.sanity(rigid, [({"plot_no": "1"}, notched, "conformed"), ({"plot_no": "2"}, hairpin, "conformed")])
    assert "2 plot(s) rigid-clipped" in notes
    assert all(p[2] == "conformed,rigid-clipped" for p in parts)
    assert parts[0][1].equals(_sq(0, 0, 20, 20)) and parts[1][1].equals(_sq(20, 0, 30, 20))
    assert not topology._spiky(parts[0][1]) and not topology._spiky(parts[1][1])


def test_sanity_removes_a_sibling_overlap_by_falling_back_the_smaller_plot():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20)), ({"plot_no": "2"}, _sq(20, 0, 24, 20))]
    parts, notes = topology.sanity(rigid, [({"plot_no": "1"}, _sq(0, 0, 20, 20), ""), ({"plot_no": "2"}, _sq(18, 0, 24, 20), "conformed")])
    assert parts[1][2] == "conformed,rigid-clipped"
    assert parts[0][1].intersection(parts[1][1]).area < 0.5


def test_sanity_is_quiet_on_clean_plots():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20))]
    parts, notes = topology.sanity(rigid, [({"plot_no": "1"}, _sq(0, 0, 20, 19.5), "conformed")])
    assert notes == [] and parts[0][2] == "conformed"


def test_sanity_keeps_the_sheet_rigid_when_its_plots_came_apart():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20)), ({"plot_no": "2"}, _sq(20, 0, 30, 20))]
    apart = [({"plot_no": "1"}, _sq(0, 0, 20, 20), "kept (conform refused)"), ({"plot_no": "2"}, _sq(21, 0, 31, 20), "conformed")]
    parts, notes = topology.sanity(rigid, apart)
    assert "plots came apart, survey kept rigid" in notes
    assert parts[1][1].equals(_sq(20, 0, 30, 20)) and parts[1][2] == "rigid (plots came apart)"


def test_sanity_trims_a_fallback_plot_against_a_filled_sibling():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 20)), ({"plot_no": "2"}, _sq(20, 0, 30, 20))]
    needle = Polygon([(20, 10), (5, 10.1), (20, 10.2)])
    spiky1 = _sq(0, 0, 20, 20).difference(needle)
    grown2 = _sq(18, 0, 30, 20)                       # sibling took a filled strip, not spiky
    parts, notes = topology.sanity(rigid, [({"plot_no": "1"}, spiky1, ""), ({"plot_no": "2"}, grown2, "filled")])
    assert parts[0][1].intersection(parts[1][1]).area < 0.5
    assert parts[0][2] == "rigid-clipped"


def test_conform_as_body_closes_small_gaps_to_settled_parcels_and_keeps_plots_clean():
    # a 20 x 10 sheet with two plots, placed 1.2 m short of a settled parcel on its east side
    rigid = [({"plot_no": "1"}, _sq(0, 0, 12, 10)), ({"plot_no": "2"}, _sq(12, 0, 20, 10))]
    settled = {"E": _sq(21.2, -5, 40, 15), "N": _sq(-5, 10.5, 25, 20)}
    parts, note = topology.conform_as_body("x", rigid, settled, anchors=set(settled), anchor_tol=1.5)
    assert parts is not None and note.startswith("conformed as one body")
    body = parts[0][1].union(parts[1][1])
    assert body.distance(settled["E"]) < 0.01 and body.distance(settled["N"]) < 0.01
    assert all(body.intersection(g).area < 0.5 for g in settled.values())
    assert parts[1][1].area > 80 and parts[1][2].startswith("body-conformed")


def test_conform_as_body_refuses_when_it_would_overlap():
    rigid = [({"plot_no": "1"}, _sq(0, 0, 20, 10))]
    settled = {"E": _sq(19.0, -5, 40, 15)}       # already overlapping by a metre
    parts, why = topology.conform_as_body("x", rigid, settled, anchors=set(settled), anchor_tol=1.5)
    assert parts is None or all(p[1].intersection(settled["E"]).area < 0.5 for p in parts)
