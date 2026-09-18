from autogeoref import georef_village


def test_parser_defaults():
    args = georef_village.parse(["35_04_077"])
    assert args.village == "35_04_077"
    assert args.raster is False and args.topology is False and args.review is True
    assert args.tracker is False, "the tracker is not written unless asked (user rule)"


def test_parser_flags():
    args = georef_village.parse(["35_04_052", "--raster", "--topology", "--no-review", "--tracker"])
    assert args.raster and args.topology and args.review is False and args.tracker


def test_refit_requires_a_survey_list():
    args = georef_village.parse(["35_04_077", "--refit", "46B,47A"])
    assert args.refit == "46B,47A"


def test_main_returns_nonzero_for_an_unknown_village(capsys):
    assert georef_village.main(["99_99_999"]) != 0
    assert "no sheets" in capsys.readouterr().out.lower()
