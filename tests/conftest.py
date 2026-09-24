"""Tests never touch the real project: a test that needs a village gets a copy of it.

On 2026-09-24 two engine tests rewrote anchors.csv and georef_status.csv in the real
Kizhikaranai folder while Akash was working on it.
"""
import shutil

import pytest

from autogeoref import paths

REAL_PROJECT = paths.PROJECT


@pytest.fixture
def village_copy(tmp_path, monkeypatch):
    """Copy FMB_Vector/<village> (and the worklist) into tmp and point paths.PROJECT there."""
    def make(village):
        src = REAL_PROJECT / "FMB_Vector" / village
        if not src.exists():
            pytest.skip("village not present")
        dst = tmp_path / "FMB_Vector" / village
        shutil.copytree(src, dst)
        wl = REAL_PROJECT / "FMB_Vector" / "worklist.csv"
        if wl.exists():
            shutil.copy2(wl, tmp_path / "FMB_Vector" / "worklist.csv")
        (tmp_path / "_logs").mkdir(exist_ok=True)
        monkeypatch.setattr(paths, "PROJECT", tmp_path)
        monkeypatch.setattr(paths, "QGIS_PROJECT", tmp_path / "no_project.qgz")
        return dst
    return make


@pytest.fixture(autouse=True)
def _real_project_is_read_only(monkeypatch):
    """Any test that forgets the copy and writes a parcel or record file into the real project fails."""
    import time
    watched = [p for p in (REAL_PROJECT / "FMB_Vector").glob("35_04_07*/*") if p.is_file()] \
        if (REAL_PROJECT / "FMB_Vector").exists() else []
    before = {p: p.stat().st_mtime for p in watched}
    yield
    changed = [p.name for p, m in before.items() if p.exists() and p.stat().st_mtime != m]
    assert not changed, "a test wrote into the real project: %s" % changed[:5]
