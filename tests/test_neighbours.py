import pytest

from autogeoref import neighbours

VILLAGE = "35_04_077"


def test_normalise_strips_prefixes_and_slashes():
    assert neighbours.normalise("103/5A") == "103"
    assert neighbours.normalise("S.No 47B") == "47B"
    assert neighbours.normalise(" 170 ") == "170"
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
