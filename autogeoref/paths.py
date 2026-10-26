"""Every filesystem path the tool uses, in one place."""
import re
from pathlib import Path

PROJECT = Path(r"D:\Projects\Tambaram_Chengalpattu")
CODE = Path(__file__).resolve().parent.parent
PUVI = Path(r"D:\Data\Puvi_data\vector")
BUFFER_GPKG = PROJECT / "Railway_Buffer_Vector_Plots.gpkg"
BUFFER_LAYER = "vector_in_buffer_30m"
RAIL_GPKG = PROJECT / "Tambaram_Chengalpattu_Railway.gpkg"
RAIL_LAYER = "rail_line"
TRACKER = PROJECT / "Georeferencing_Tracker.xlsx"
# Akash's QGIS project: the layers it points at are the parcels as he sees them (see visible.py)
QGIS_PROJECT = Path(r"D:\Projects\CUMTA\15Sep26\New folder\CUMTA_Georeferencing.qgz")


def vector_dir(village):
    return PROJECT / "FMB_Vector" / village


def sketch_dir(village):
    return PROJECT / "FMB_Sketches" / village


def georef_dir(village):
    return PROJECT / "FMB_Georef" / village


def logs_dir():
    return PROJECT / "_logs"


def sheet_path(village, survey):
    return vector_dir(village) / ("%s_parcels.geojson" % survey)


def points_path(village, survey):
    return vector_dir(village) / ("%s_parcels.geojson.points" % survey)


def auto_dir(village):
    """The tool's own folder. Nothing in it is Akash's; nothing of his is ever written from here."""
    return vector_dir(village) / "_auto"


def output_path(village, survey):
    """Where the engine writes a placement: its own _auto folder, never a hand file name.

    Until 2026-09-24 this was <s>_parcels_modified.gpkg in his folder, so tool output posed as
    hand work and approved parcels ended up in two files (see visible.py).
    """
    return auto_dir(village) / ("%s_parcels_auto.gpkg" % survey)


def status_path(village):
    return vector_dir(village) / "georef_status.csv"


def anchors_path(village):
    return vector_dir(village) / "anchors.csv"


def gcp_csv(village, survey):
    return vector_dir(village) / "gcp" / ("%s_gcp.csv" % survey)


def raster_pin(village):
    return georef_dir(village) / "raster.json"


def puvi_shapefile(village):
    """D:\\Data\\Puvi_data\\vector\\35\\04\\077\\vector\\*.shp for village 35_04_077."""
    d, t, v = village.split("_")
    hits = sorted((PUVI / d / t / v / "vector").glob("*.shp"))
    return hits[0] if hits else None


def survey_sort_key(survey):
    """'40A' -> (40, '40A') so 9 < 40A < 40B < 169."""
    m = re.match(r"\d+", str(survey))
    return (int(m.group()) if m else 10 ** 9, str(survey))


def manual_files(village, survey):
    """Every hand-placed GeoPackage for this survey, the one his QGIS project shows first.

    Newest modification time is only the fallback: on 2026-09-23 an approved parcel lived in
    <s>_parcels_auto.gpkg (shown in QGIS) and a stale copy in <s>_parcels_modified.gpkg (edited by
    the tool), and "newest" picked the copy he never saw.
    """
    hits = sorted(vector_dir(village).glob("%s_parcels_modified*.gpkg" % survey),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    try:
        from . import visible
        seen = visible.shown(village).get(str(survey))
    except Exception:
        seen = None
    if seen is not None:
        hits = [seen[0]] + [h for h in hits if h.resolve() != seen[0]]
    return hits


def surveys_with_sheets(village):
    """Survey ids that have a sheet GeoJSON, numerically sorted. An empty file (a conversion
    that drew nothing, like 52A on 2026-09-23) is not a sheet."""
    out = [p.name[: -len("_parcels.geojson")] for p in vector_dir(village).glob("*_parcels.geojson")
           if p.stat().st_size > 200]
    return sorted(out, key=survey_sort_key)
