"""Which file, and which layer inside it, QGIS actually shows for a survey. The one source of truth.

Akash's audit of 2026-09-24 found that since 23 Sep the tool and his QGIS had been working on
different copies of the same parcels:

  * his hand GeoPackages keep their parcels in a layer named after the file
    (565_parcels_modified.gpkg -> layer 565_parcels_modified); scripts wrote layer "parcels" into
    the same file, which ADDED a second, hidden layer instead of changing his, and every later
    edit, check and export ran on that hidden copy;
  * tool-written files carry an "edges" line layer, and in 21 files it came first, so a
    read_file() without a layer name loaded boundary lines as if they were the parcel;
  * approved parcels exist twice: QGIS shows <s>_parcels_auto.gpkg, the tool edited
    <s>_parcels_modified.gpkg.

So nothing reads a parcel GeoPackage without naming the layer, the layer is chosen by
`hand_layer`, and the file for a survey is the one the QGIS project points at (`shown`).
"""
import functools
import re
import zipfile
from pathlib import Path

import geopandas as gpd
import pyogrio

from . import paths

TOOL_LAYER = "parcels"
LINE_LAYER = "edges"


def polygon_layers(path):
    """Names of the polygon layers in a GeoPackage, in file order."""
    return [str(n) for n, g in pyogrio.list_layers(str(path))
            if g is not None and "polygon" in str(g).lower()]


def hand_layer(path):
    """The layer QGIS shows in a parcel GeoPackage.

    His own layer carries the file's name; a file the tool wrote has only "parcels". A file with
    both is his file with a hidden tool copy in it: his layer wins. Never the "edges" lines.
    """
    names = polygon_layers(path)
    stem = Path(path).stem
    if stem in names:
        return stem
    if TOOL_LAYER in names:
        return TOOL_LAYER
    if len(names) == 1:
        return names[0]
    raise ValueError("%s: cannot tell which polygon layer QGIS shows (%s)" % (path, names))


def read_hand(path, **kw):
    """The parcels exactly as QGIS shows them."""
    return gpd.read_file(str(path), layer=hand_layer(path), **kw)


def has_hidden_copy(path):
    """True when a hand file also holds a tool-written "parcels" layer QGIS does not show."""
    names = polygon_layers(path)
    return Path(path).stem in names and TOOL_LAYER in names


def refuse_tool_write(path):
    """Raise before the tool writes its "parcels" layer into a file where QGIS shows another layer.

    Into his hand file that write either ADDS a hidden copy beside his layer (what happened in 27
    files on 23-24 Sep) or, through a whole-file replace, deletes his layer. Neither is allowed.
    """
    p = Path(path)
    if is_shown(p):
        raise RuntimeError("refusing to write into %s: his QGIS project shows it" % p.name)
    if p.exists() and polygon_layers(p) and hand_layer(p) != TOOL_LAYER:
        raise RuntimeError("refusing to write layer %r into %s: QGIS shows %r there (a hand file)"
                           % (TOOL_LAYER, p.name, hand_layer(p)))


def is_shown(path, qgz=None):
    """True when his QGIS project shows this file (in any village): it is his, the tool keeps out."""
    try:
        target = Path(path).resolve()
        return any(fp == target for _n, fp, _l in project_sources(qgz))
    except (OSError, IndexError, KeyError, zipfile.BadZipFile):
        return False


def project_sources(qgz=None):
    """[(layer name in QGIS, absolute file path, layername)] for every GeoPackage layer of the project."""
    qgz = Path(qgz or paths.QGIS_PROJECT)
    return list(_project_sources(str(qgz), qgz.stat().st_mtime))


@functools.lru_cache(maxsize=4)
def _project_sources(qgz, _mtime):
    qgz = Path(qgz)
    with zipfile.ZipFile(qgz) as z:
        xml = z.read([n for n in z.namelist() if n.endswith(".qgs")][0]).decode("utf-8")
    out = []
    for block in re.findall(r"<maplayer\b.*?</maplayer>", xml, re.S):
        ds = re.search(r"<datasource>(.*?)</datasource>", block, re.S)
        nm = re.search(r"<layername>(.*?)</layername>", block, re.S)
        if not ds or ".gpkg" not in ds.group(1):
            continue
        src = ds.group(1).replace("&amp;", "&")
        f, _sep, rest = src.partition("|")
        m = re.search(r"layername=([^|]+)", rest)
        fp = Path(f)
        if not fp.is_absolute():
            fp = (qgz.parent / fp).resolve()
        out.append((nm.group(1) if nm else "", fp, m.group(1) if m else None))
    return tuple(out)


def shown(village, qgz=None):
    """{survey: (file, layer)} for the parcel GeoPackages of `village` that the QGIS project shows."""
    vd = paths.vector_dir(village).resolve()
    out = {}
    for _name, fp, layer in project_sources(qgz):
        if fp.parent != vd:
            continue
        m = re.match(r"(.+?)_parcels_", fp.name)
        if not m:
            continue
        s = m.group(1)
        if s in out and out[s] != (fp, layer):
            raise ValueError("%s %s: the project shows two files (%s, %s)" % (village, s, out[s][0].name, fp.name))
        out[s] = (fp, layer or hand_layer(fp))
    return out
