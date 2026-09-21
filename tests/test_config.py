import json

import pytest

from autogeoref import config


def test_defaults_are_the_rule_based_baseline():
    cfg = config.load()
    assert cfg.ml.confidence == "rule"
    assert cfg.imagery.source == "google_xyz"
    assert cfg.acceptance_m == 3.0
    assert cfg.stretch.rule in config.STRETCH_RULES


def test_unknown_key_is_an_error_not_a_silent_default(tmp_path):
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"acceptance_m": 3.0, "imagery": {"sourse": "google_xyz"}}), encoding="utf-8")
    with pytest.raises(ValueError, match="sourse"):
        config.load(str(p))


def test_environment_variable_overrides_the_default(tmp_path, monkeypatch):
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"acceptance_m": 2.5, "stretch": {"rule": "pct", "pct": 7.5}}), encoding="utf-8")
    monkeypatch.setenv(config.ENV_VAR, str(p))
    cfg = config.load()
    assert cfg.acceptance_m == 2.5 and cfg.stretch.rule == "pct" and cfg.stretch.pct == 7.5
    assert cfg.stretch.band_pct == 5.0, "unspecified keys keep their defaults"
    assert cfg.source_path == str(p)


def test_invalid_switches_are_rejected(tmp_path):
    p = tmp_path / "c.json"
    p.write_text(json.dumps({"ml": {"confidence": "unet"}}), encoding="utf-8")
    with pytest.raises(ValueError):
        config.load(str(p))
    p.write_text(json.dumps({"imagery": {"source": "bing"}}), encoding="utf-8")
    with pytest.raises(ValueError):
        config.load(str(p))


def test_cache_dir_is_outside_the_repository_and_under_the_project(monkeypatch, tmp_path):
    from autogeoref import paths
    monkeypatch.setattr(paths, "PROJECT", tmp_path)
    cfg = config.Config()
    d = config.cache_dir(cfg)
    assert d == tmp_path / "_cache" / "imagery"
    assert str(paths.CODE) not in str(d)
    cfg.imagery.cache_dir = str(tmp_path / "elsewhere")
    assert config.cache_dir(cfg) == tmp_path / "elsewhere"
