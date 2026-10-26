from pathlib import Path

from autogeoref import paths


def test_village_directories_point_at_the_project():
    assert paths.vector_dir("35_04_077") == paths.PROJECT / "FMB_Vector" / "35_04_077"
    assert paths.sketch_dir("35_04_077") == paths.PROJECT / "FMB_Sketches" / "35_04_077"
    assert paths.georef_dir("35_04_077") == paths.PROJECT / "FMB_Georef" / "35_04_077"


def test_sheet_and_output_paths():
    assert paths.sheet_path("35_04_077", "47B").name == "47B_parcels.geojson"
    assert paths.points_path("35_04_077", "47B").name == "47B_parcels.geojson.points"
    assert paths.output_path("35_04_077", "47B").name == "47B_parcels_auto.gpkg"
    assert paths.output_path("35_04_077", "47B").parent.name == "_auto"   # the tool's own folder


def test_manual_files_put_the_file_qgis_shows_first():
    files = paths.manual_files("35_04_077", "47A")
    assert [f.name for f in files][:1] == ["47A_parcels_modified2.gpkg"]
    # 24 Sep clean-up: one file per survey; older variants live in _recovery_20260924/archive
    assert len(files) == 1


def test_surveys_sort_numerically_then_by_letter():
    got = sorted(["169", "40B", "42A", "9", "40A"], key=paths.survey_sort_key)
    assert got == ["9", "40A", "40B", "42A", "169"]


def test_surveys_with_sheets_finds_the_village():
    got = paths.surveys_with_sheets("35_04_077")
    assert "47B" in got and "169" in got
    assert all(paths.sheet_path("35_04_077", s).exists() for s in got)
