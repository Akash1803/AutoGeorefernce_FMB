"""Export the project's Google Satellite XYZ layer to a georeferenced GeoTIFF and pin it.

The QGIS layer "Google Satellite" is the XYZ source
https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z} in EPSG:3857. rasterio reads it through GDAL's
WMS/TMS driver given the XML description below, so no screenshot and no ECW writer are needed.
"""
import datetime
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import Resampling, reproject, transform_bounds

from . import paths

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
</GDAL_WMS>
"""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def export(village, bounds, zoom=20, resolution=0.15, run_id=None):
    """Write the satellite window for `bounds` (EPSG:32644 minx, miny, maxx, maxy) as a GeoTIFF."""
    out_dir = paths.georef_dir(village)
    out_dir.mkdir(parents=True, exist_ok=True)
    desc = out_dir / ("google_sat_z%d.xml" % zoom)
    desc.write_text(TILE_XML.replace("<TileLevel>20</TileLevel>", "<TileLevel>%d</TileLevel>" % zoom),
                    encoding="utf-8")
    minx, miny, maxx, maxy = bounds
    width = int(round((maxx - minx) / resolution))
    height = int(round((maxy - miny) / resolution))
    dst_tr = from_origin(minx, maxy, resolution, resolution)
    dst = np.zeros((3, height, width), np.uint8)
    with rasterio.open(str(desc)) as src:
        win = src.window(*transform_bounds("EPSG:32644", src.crs, minx, miny, maxx, maxy))
        win = win.round_offsets().round_lengths()
        arr = src.read(window=win)
        src_tr = src.window_transform(win)
        for band in range(3):
            reproject(arr[band], dst[band], src_transform=src_tr, src_crs=src.crs,
                      dst_transform=dst_tr, dst_crs="EPSG:32644", resampling=Resampling.bilinear)
    stamp = run_id or datetime.datetime.now().strftime("%Y%m%d")
    out = out_dir / ("satellite_z%d_%s.tif" % (zoom, stamp))
    previous = pinned(village)
    if previous and previous.get("path") and str(out) != previous["path"]:
        archive = paths.logs_dir() / ("georef_raster_archive_%s" % stamp)
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
    paths.raster_pin(village).write_text(json.dumps(
        {"path": str(out), "sha256": _sha256(out), "zoom": zoom, "resolution_m": resolution,
         "fetched": datetime.datetime.now().isoformat(timespec="seconds"),
         "bounds_32644": [minx, miny, maxx, maxy],
         "source": "QGIS layer 'Google Satellite' (XYZ mt1.google.com/vt/lyrs=s) via GDAL_WMS/TMS",
         "note": "internal review basemap only; not for redistribution"}, indent=1), encoding="utf-8")
    return out


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
