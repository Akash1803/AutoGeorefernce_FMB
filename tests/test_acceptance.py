import pytest

from autogeoref import anchors, engine, paths, review

pytestmark = pytest.mark.slow
VILLAGE = "35_04_077"


def _skip_if_absent():
    if not paths.vector_dir(VILLAGE).exists():
        pytest.skip("village not present")


def test_a_full_run_leaves_every_manual_file_and_points_file_byte_identical():
    _skip_if_absent()
    manual = {p: p.read_bytes() for p in paths.vector_dir(VILLAGE).glob("*_parcels_modified*.gpkg")}
    points = {p: p.read_bytes() for p in paths.vector_dir(VILLAGE).glob("*.points")}
    engine.run(VILLAGE, do_raster=False, do_topology=True, do_review=False)
    for p, data in manual.items():
        assert p.read_bytes() == data, "%s was modified" % p.name
    for p, data in points.items():
        assert p.read_bytes() == data, "%s was modified" % p.name


def test_every_buffer_parcel_is_classified():
    _skip_if_absent()
    engine.run(VILLAGE, do_raster=False, do_topology=False, do_review=False)
    rows = review.read_status(VILLAGE)
    got = {r["survey"] for r in rows} | set(anchors.load_anchors(VILLAGE, set()))
    assert set(paths.surveys_with_sheets(VILLAGE)) <= got
    assert all(r.get("colour") in ("green", "amber", "red") for r in rows)


def test_railway_strips_are_not_flipped_when_they_are_placed():
    """169, 170 and 171 are near-symmetric; only the transcription can orient them."""
    _skip_if_absent()
    rows = {r["survey"]: r for r in review.read_status(VILLAGE)}
    for s in ("169", "170", "171"):
        if s in rows and rows[s].get("status") == "placed":
            assert rows[s]["colour"] != "green" or rows[s].get("observable") == "True"
