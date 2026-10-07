"""Satellite check of placed (or proposed) parcels: how much of the outline runs along a straight line
visible in the Google Satellite imagery, and whether the outline would fit better shifted.

It never moves a parcel - on 2026-09-20 the imagery preferred a road edge over the true boundary - it only
reports, next to the same measure on Akash's own parcels in the same window, so poor imagery is not
blamed on the parcel.
"""
import math
import os

import numpy as np

from . import align, edges, raster, sheets

ZOOM, RES_M, MARGIN_M = 19, 0.3, 40.0
SHIFT_M, STEP_M = 6.0, 0.5
FLAG_SHIFT_M = 2.0          # a better fit this far away is reported
FLAG_GAIN = 0.05            # ... when it gains at least this share of the outline


MAX_PX = 4000               # longest side of the rendered window; resolution coarsens beyond it


def render_from_qgis(bounds, out_png, layer_name="Google Satellite"):
    """Render the project's satellite layer (QGIS has the tiles cached) to a PNG + world file in
    EPSG:32644. A 100 m window took 24 s through GDAL's tile reader on 2026-10-07; QGIS takes seconds."""
    from . import qgis_bridge
    minx, miny, maxx, maxy = (float(v) for v in bounds)       # numpy floats repr as np.float64(...)
    res = max(RES_M, max(maxx - minx, maxy - miny) / MAX_PX)
    w, h = int(round((maxx - minx) / res)), int(round((maxy - miny) / res))
    code = """
import json
from qgis.core import QgsProject, QgsMapSettings, QgsMapRendererParallelJob, QgsRectangle, QgsCoordinateReferenceSystem
from qgis.PyQt.QtCore import QSize
from qgis.PyQt.QtGui import QColor
L = QgsProject.instance().mapLayersByName(%r)[0]
ms = QgsMapSettings(); ms.setLayers([L]); ms.setDestinationCrs(QgsCoordinateReferenceSystem("EPSG:32644"))
ms.setExtent(QgsRectangle(%r, %r, %r, %r)); ms.setOutputSize(QSize(%d, %d)); ms.setBackgroundColor(QColor(0, 0, 0))
j = QgsMapRendererParallelJob(ms); j.start(); j.waitForFinished()
ok = j.renderedImage().save(%r, "JPG", 95)
print(json.dumps({"ok": ok}))
""" % (layer_name, minx, miny, maxx, maxy, w, h, str(out_png).replace(os.sep, "/"))
    if not qgis_bridge.json_result(code, timeout=300).get("ok"):
        raise RuntimeError("QGIS could not render the satellite window")
    with open(str(out_png)[:-4] + ".jgw", "w") as fh:
        fh.write("\n".join(["%.6f" % res, "0", "0", "%.6f" % -res, "%.6f" % (minx + res / 2), "%.6f" % (maxy - res / 2)]) + "\n")
    return out_png


class Imagery:
    def __init__(self, village, polys, run_id, out_dir=None):
        b = np.array([p.bounds for p in polys])
        bounds = (b[:, 0].min() - MARGIN_M, b[:, 1].min() - MARGIN_M, b[:, 2].max() + MARGIN_M, b[:, 3].max() + MARGIN_M)
        if out_dir is not None:
            self.tif = render_from_qgis(bounds, out_dir / "satellite_window.jpg")
        else:
            self.tif = raster.export(village, bounds, zoom=ZOOM, resolution=RES_M, run_id=run_id)
        self.index = edges.SegmentIndex(edges.detect(self.tif, min_len_m=3.0))
        self.bounds = bounds

    def _samples(self, poly):
        ring = list(poly.exterior.coords)[:-1]
        return sheets.sample_outline(ring, 1.0)

    def share(self, poly, dx=0.0, dy=0.0):
        pts, brg = self._samples(poly)
        return align.share(pts, brg, self.index, 0.0, (dx, dy))

    def shift_test(self, poly):
        """(share here, best share within +-SHIFT_M, its offset dx, dy)."""
        pts, brg = self._samples(poly)
        here = align.share(pts, brg, self.index, 0.0, (0.0, 0.0))
        best = (here, 0.0, 0.0)
        for dx in np.arange(-SHIFT_M, SHIFT_M + 1e-9, STEP_M):
            for dy in np.arange(-SHIFT_M, SHIFT_M + 1e-9, STEP_M):
                s = align.share(pts, brg, self.index, 0.0, (dx, dy))
                if s > best[0] + 1e-9:
                    best = (s, float(dx), float(dy))
        return here, best[0], best[1], best[2]

    def verdict(self, poly, baseline):
        here, best, dx, dy = self.shift_test(poly)
        d = math.hypot(dx, dy)
        out = {"sat_share": round(here, 2), "sat_share_his_parcels": round(baseline, 2),
               "sat_best_shift_m": round(d, 1), "sat_best_share": round(best, 2)}
        if d >= FLAG_SHIFT_M and best - here >= FLAG_GAIN:
            bearing = (math.degrees(math.atan2(dx, dy)) + 360) % 360
            out["sat_note"] = "imagery fits better %.1f m towards %s" % (d, _compass(bearing))
            out["sat_flag"] = "red" if d >= 4.0 and best - here >= 2 * FLAG_GAIN else "amber"
        else:
            out["sat_note"] = "fits the imagery where it is"
            out["sat_flag"] = "green" if here >= 0.5 * baseline else "amber"
            if out["sat_flag"] == "amber":
                out["sat_note"] = "few visible lines along its boundary (%.2f vs %.2f on your parcels)" % (here, baseline)
        return out


def _compass(b):
    return ["N", "NE", "E", "SE", "S", "SW", "W", "NW"][int(((b + 22.5) % 360) // 45)]
