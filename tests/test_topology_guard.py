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
