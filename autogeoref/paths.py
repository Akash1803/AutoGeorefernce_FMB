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


def output_path(village, survey):
    return vector_dir(village) / ("%s_parcels_modified.gpkg" % survey)


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
    """Every hand-placed GeoPackage for this survey, newest modification time first."""
    hits = list(vector_dir(village).glob("%s_parcels_modified*.gpkg" % survey))
    return sorted(hits, key=lambda p: p.stat().st_mtime, reverse=True)


def surveys_with_sheets(village):
    """Survey ids that have a sheet GeoJSON, numerically sorted. An empty file (a conversion
    that drew nothing, like 52A on 2026-09-23) is not a sheet."""
    out = [p.name[: -len("_parcels.geojson")] for p in vector_dir(village).glob("*_parcels.geojson")
           if p.stat().st_size > 200]
    return sorted(out, key=survey_sort_key)
