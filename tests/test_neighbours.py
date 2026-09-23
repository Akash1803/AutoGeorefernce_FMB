import pytest

from autogeoref import neighbours

VILLAGE = "35_04_077"


def test_normalise_strips_prefixes_and_slashes():
    assert neighbours.normalise("103/5A") == "103"
    assert neighbours.normalise("S.No 47B") == "47B"
    assert neighbours.normalise(" 170 ") == "170"
    assert neighbours.normalise("569/C") == "569C", "a slash letter is the portal's subdivision suffix"
    assert neighbours.normalise("569/2") == "569", "a slash number is a plot inside the survey"
    assert neighbours.normalise("V.No. 74 THIRUKACHUR") is None, "a village number is not a survey"
    assert neighbours.normalise("Railway") is None


def test_printed_numbers_match_the_transcription():
    got = neighbours.printed(VILLAGE, "171")
    assert {"170", "48A", "48B", "47B", "47A", "46A", "46B"} <= got


def test_mutual_naming_makes_a_neighbour_pair():
    assert neighbours.are_neighbours(VILLAGE, "171", "170") is True
    assert neighbours.are_neighbours(VILLAGE, "43A", "46B") is True


def test_two_fully_read_sheets_that_ignore_each_other_are_not_neighbours():
    assert neighbours.are_neighbours(VILLAGE, "40A", "171") is False


def test_unknown_when_a_village_has_no_transcription():
    assert neighbours.are_neighbours("35_04_052", "69", "68") is None


def test_side_and_side_vector():
    assert neighbours.side(VILLAGE, "171", "170") == "S"
    assert neighbours.side_vector("N") == (0.0, 1.0)
    assert neighbours.side_vector("SE") == pytest.approx((0.7071, -0.7071), abs=1e-3)



def test_a_read_sheet_is_the_complete_list_of_its_neighbours(tmp_path, monkeypatch):
    import json
    from autogeoref import paths
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    d = tmp_path / "FMB_Vector" / "V" / "neighbour_transcription"; d.mkdir(parents=True)
    (d / "nb_override_V.json").write_text(json.dumps({"613": [{"number": "565", "side": "E", "readers": 2}]}), encoding="utf-8")
    neighbours._CACHE.pop("V", None)
    assert neighbours.are_neighbours("V", "613", "565") is True
    assert neighbours.are_neighbours("V", "613", "609") is False, "613 was read and does not print 609"
    assert neighbours.are_neighbours("V", "609", "613") is False
    assert neighbours.are_neighbours("V", "609", "610") is None, "neither sheet was read"
