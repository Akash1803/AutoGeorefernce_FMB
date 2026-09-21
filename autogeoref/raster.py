"""Export the project's Google Satellite XYZ layer to a georeferenced GeoTIFF and pin it.

The QGIS layer "Google Satellite" is the XYZ source
https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z} in EPSG:3857. rasterio reads it through GDAL's
WMS/TMS driver given the XML description below, so no screenshot and no ECW writer are needed.
"""
import csv
import datetime
import time
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds

from . import config as configmod, paths

LOG_COLUMNS = ["fetched", "village", "source", "zoom", "resolution_m", "minx", "miny", "maxx", "maxy",
               "sha256", "path", "run_id"]

TILE_XML = """<GDAL_WMS>
  <Service name="TMS"><ServerUrl>https://mt1.google.com/vt/lyrs=s&amp;x=${x}&amp;y=${y}&amp;z=${z}</ServerUrl></Service>
  <DataWindow>
    <UpperLeftX>-20037508.34</UpperLeftX><UpperLeftY>20037508.34</UpperLeftY>
    <LowerRightX>20037508.34</LowerRightX><LowerRightY>-20037508.34</LowerRightY>
    <TileLevel>20</TileLevel><TileCountX>1</TileCountX><TileCountY>1</TileCountY><YOrigin>top</YOrigin>
  </DataWindow>
  <Projection>EPSG:3857</Projection>
  <BlockSizeX>256</BlockSizeX><BlockSizeY>256</BlockSizeY><BandsCount>3</BandsCount>
  <MaxConnections>4</MaxConnections>
  <Timeout>60</Timeout>
  <ZeroBlockHttpCodes>204,403,404,429,500,502,503,504</ZeroBlockHttpCodes>
  <ZeroBlockOnServerException>true</ZeroBlockOnServerException>
</GDAL_WMS>
"""
READ_ATTEMPTS = 3          # one refused tile must not abort a 300 m window (2026-09-20)
MAX_EMPTY_SHARE = 0.02     # more empty pixels than this and the export is refused, not pinned


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def export(village, bounds, zoom=20, resolution=0.15, run_id=None, cfg=None):
    """Write the imagery window for `bounds` (EPSG:32644 minx, miny, maxx, maxy) as a GeoTIFF.

    The window goes to the imagery cache (gitignored, deletable, see ``config.cache_dir``), never
    under the repository, and every export is appended to ``imagery_log.csv`` there with its
    source. ``cfg.imagery.source`` is ``google_xyz`` (tiles, inference only) or ``geotiff:<path>``
    (a licensed orthomosaic); the pin file records which one produced the window.
    """
    cfg = cfg or configmod.load()
    zoom = zoom if zoom is not None else cfg.imagery.zoom
    out_dir = configmod.cache_dir(cfg) / village
    out_dir.mkdir(parents=True, exist_ok=True)
    paths.georef_dir(village).mkdir(parents=True, exist_ok=True)
    source = cfg.imagery.source
    if source == "google_xyz":
        desc = out_dir / ("google_sat_z%d.xml" % zoom)
        desc.write_text(TILE_XML.replace("<TileLevel>20</TileLevel>", "<TileLevel>%d</TileLevel>" % zoom),
                        encoding="utf-8")
    else:
        desc = Path(source[len("geotiff:"):])
        if not desc.exists():
            raise FileNotFoundError("imagery source %s not found" % desc)
    minx, miny, maxx, maxy = bounds
    width = int(round((maxx - minx) / resolution))
    height = int(round((maxy - miny) / resolution))
    dst_tr = from_origin(minx, maxy, resolution, resolution)
    dst = np.zeros((3, height, width), np.uint8)
    last = None
    for attempt in range(READ_ATTEMPTS):
        try:
            with rasterio.open(str(desc)) as src:
                win = src.window(*transform_bounds("EPSG:32644", src.crs, minx, miny, maxx, maxy))
                win = win.round_offsets().round_lengths()
                arr = src.read(window=win)
                src_tr = src.window_transform(win)
                src_crs = src.crs
            empty = float((arr.max(axis=0) == 0).mean())
            if empty <= MAX_EMPTY_SHARE:
                break
            last = RuntimeError("%.1f %% of the window came back empty" % (100 * empty))
        except Exception as exc:                      # a refused tile, a dropped connection
            last = exc
        time.sleep(3.0 * (attempt + 1))
    else:
        raise RuntimeError("satellite export failed after %d attempts: %s" % (READ_ATTEMPTS, last))
    for band in range(3):
        reproject(arr[band], dst[band], src_transform=src_tr, src_crs=src_crs,
                  dst_transform=dst_tr, dst_crs="EPSG:32644", resampling=Resampling.bilinear)
    if float((dst.max(axis=0) == 0).mean()) > MAX_EMPTY_SHARE:
        raise RuntimeError("satellite window is empty after resampling; nothing pinned")
    stamp = run_id or datetime.datetime.now().strftime("%Y%m%d")
    out = out_dir / ("satellite_z%d_%s.tif" % (zoom, stamp))
    previous = pinned(village)
    if previous and previous.get("path") and str(out) != previous["path"]:
        archive = out_dir / "_archive"
        archive.mkdir(parents=True, exist_ok=True)
        try:
            shutil.move(previous["path"], str(archive / Path(previous["path"]).name))
        except (OSError, FileNotFoundError):
            pass
    profile = {"driver": "GTiff", "height": height, "width": width, "count": 3, "dtype": "uint8",
               "crs": "EPSG:32644", "transform": dst_tr, "compress": "DEFLATE", "tiled": True}
    with rasterio.open(out, "w", **profile) as fh:
        fh.write(dst)
        fh.build_overviews([2, 4, 8, 16], Resampling.average)
    sha = _sha256(out)
    fetched = datetime.datetime.now().isoformat(timespec="seconds")
    paths.raster_pin(village).write_text(json.dumps(
        {"path": str(out), "sha256": sha, "zoom": zoom, "resolution_m": resolution,
         "fetched": fetched, "bounds_32644": [minx, miny, maxx, maxy],
         "source": source, "cache_dir": str(configmod.cache_dir(cfg)),
         "note": "inference only; internal review; not for redistribution; delete with --purge-cache"},
        indent=1), encoding="utf-8")
    log_path = configmod.cache_dir(cfg) / "imagery_log.csv"
    new_log = not log_path.exists()
    with log_path.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_COLUMNS)
        if new_log:
            w.writeheader()
        w.writerow({"fetched": fetched, "village": village, "source": source, "zoom": zoom,
                    "resolution_m": resolution, "minx": minx, "miny": miny, "maxx": maxx, "maxy": maxy,
                    "sha256": sha, "path": str(out), "run_id": run_id or ""})
    return out


def purge_cache(cfg=None):
    """Delete every cached imagery window and its log; pins are left and become stale."""
    cfg = cfg or configmod.load()
    d = configmod.cache_dir(cfg)
    if d.exists():
        shutil.rmtree(d)
    return d


def migrate_to_cache(village, cfg=None):
    """Move a window pinned under FMB_Georef (pre-PR 0 layout) into the cache and repoint the pin."""
    cfg = cfg or configmod.load()
    pin = pinned(village)
    if not pin or not pin.get("path"):
        return None
    src = Path(pin["path"])
    dst_dir = configmod.cache_dir(cfg) / village
    if not src.exists() or dst_dir in src.parents:
        return src if src.exists() else None
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    shutil.move(str(src), str(dst))
    pin["path"] = str(dst)
    pin.setdefault("source", "google_xyz")
    pin["cache_dir"] = str(configmod.cache_dir(cfg))
    paths.raster_pin(village).write_text(json.dumps(pin, indent=1), encoding="utf-8")
    return dst


def pinned(village):
    p = paths.raster_pin(village)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def open_pinned(village):
    pin = pinned(village)
    if not pin:
        raise FileNotFoundError("no raster pinned for %s; run with --raster first" % village)
    if _sha256(pin["path"]) != pin["sha256"]:
        raise ValueError("pinned raster changed on disk: %s" % pin["path"])
    return rasterio.open(pin["path"])
