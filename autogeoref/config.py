"""Run configuration for the pipeline and the evaluation harness.

One JSON file, ``configs/default.json`` in the repository, overridable with the environment
variable ``AUTOGEOREF_CONFIG=<path>``. Unknown keys are an error, so a typo cannot fall back to a
default in silence. Every ML component sits behind a switch here and every switch defaults to the
rule-based baseline.
"""
import dataclasses
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Type

from . import paths

log = logging.getLogger(__name__)

DEFAULT_PATH = paths.CODE / "configs" / "default.json"
ENV_VAR = "AUTOGEOREF_CONFIG"
STRETCH_RULES = ("pct", "axis", "displacement")
CONFIDENCE_MODES = ("rule", "gbm")


@dataclass
class ImageryConfig:
    """Where satellite pixels come from and where windows may be cached.

    ``source`` is ``google_xyz`` (Google Satellite tiles, inference only, never stored as training
    data) or ``geotiff:<path>`` for a licensed orthomosaic. ``cache_dir`` empty means
    ``<project>/_cache/imagery``; the directory is gitignored and deletable.
    """
    source: str = "google_xyz"
    cache_dir: str = ""
    zoom: int = 20
    resolution_m: float = 0.15


@dataclass
class StretchConfig:
    """When a hand placement is too stretched to be a reference for a label.

    ``rule``: ``pct`` (area scale off 1 by more than ``pct`` percent), ``axis`` (either side of
    the minimum rotated rectangle off by more than ``pct`` percent) or ``displacement`` (mean
    distance between the FMB-exact outline and the hand outline above ``displacement_mean_m``).
    Whatever the rule, a rigid fit worse than ``rms_m`` excludes. ``band_pct`` marks the middle
    tier that stays labelled but is reported separately.
    """
    rule: str = "pct"
    pct: float = 10.0
    band_pct: float = 5.0
    displacement_mean_m: float = 3.0
    rms_m: float = 3.0


@dataclass
class MLConfig:
    """Which components run: ``rule`` is the baseline colour rule; ``gbm`` arrives with PR 3."""
    confidence: str = "rule"


@dataclass
class Config:
    acceptance_m: float = 3.0
    imagery: ImageryConfig = field(default_factory=ImageryConfig)
    stretch: StretchConfig = field(default_factory=StretchConfig)
    ml: MLConfig = field(default_factory=MLConfig)
    source_path: str = ""


NESTED: Dict[str, Type] = {"imagery": ImageryConfig, "stretch": StretchConfig, "ml": MLConfig}


def _build(cls: Type, data: Dict[str, Any], where: str) -> Any:
    names = {f.name for f in dataclasses.fields(cls)}
    unknown = sorted(set(data) - names)
    if unknown:
        raise ValueError("unknown key(s) in %s: %s" % (where, ", ".join(unknown)))
    kwargs: Dict[str, Any] = {}
    for name, value in data.items():
        if name in NESTED and cls is Config:
            kwargs[name] = _build(NESTED[name], dict(value), "%s.%s" % (where, name))
        else:
            kwargs[name] = value
    return cls(**kwargs)


def validate(cfg: Config) -> Config:
    if cfg.stretch.rule not in STRETCH_RULES:
        raise ValueError("stretch.rule must be one of %s" % (STRETCH_RULES,))
    if cfg.ml.confidence not in CONFIDENCE_MODES:
        raise ValueError("ml.confidence must be one of %s" % (CONFIDENCE_MODES,))
    src = cfg.imagery.source
    if src != "google_xyz" and not src.startswith("geotiff:"):
        raise ValueError("imagery.source must be 'google_xyz' or 'geotiff:<path>'")
    if cfg.acceptance_m <= 0:
        raise ValueError("acceptance_m must be positive")
    return cfg


def load(path: Optional[str] = None) -> Config:
    """The configuration in force: explicit path, else the environment variable, else the default."""
    p = Path(path or os.environ.get(ENV_VAR) or DEFAULT_PATH)
    if not p.exists():
        log.warning("config %s not found; using built-in defaults", p)
        return validate(Config(source_path=""))
    data = json.loads(p.read_text(encoding="utf-8"))
    cfg = _build(Config, data, "config")
    cfg.source_path = str(p)
    log.info("config loaded from %s", p)
    return validate(cfg)


def cache_dir(cfg: Config) -> Path:
    """The imagery cache directory: gitignored, outside the repository, safe to delete."""
    return Path(cfg.imagery.cache_dir) if cfg.imagery.cache_dir else paths.PROJECT / "_cache" / "imagery"


def to_dict(cfg: Config) -> Dict[str, Any]:
    return dataclasses.asdict(cfg)
