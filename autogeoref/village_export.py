"""Village deliverables built from the parcels exactly as Akash's QGIS shows them.

    python -m autogeoref.village_export 35_04_074 Thirukatchur --out <folder>

Writes, in EPSG:4326:
  <Name>_Manual_Georeferenced.geojson  one row per plot: kide / survey_no / subdiv_no / area_acre
                                       first (his mandatory fields), then area_sqm, status (from
                                       his tracker), iou_vs_pdf of the survey, source_file
  fabric_check_<village>.csv           the topology/shape check of the same parcels

The 24 Sep GeoJSONs were built from hidden tool copies; this module reads only visible.shown().
"""
import argparse
import re
import zipfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

from . import fabric_check, fmb_finish, paths, sheets, visible

ACRE = 4046.8564224


def tracker_status(village, workbook=None):
    """{survey: Status} for one village, read from his tracker (read-only)."""
    workbook = Path(workbook or paths.TRACKER)
    if not workbook.exists():
        return {}
    with zipfile.ZipFile(workbook) as z:
        sheet = z.read("xl/worksheets/sheet1.xml").decode("utf-8")
        shared = z.read("xl/sharedStrings.xml").decode("utf-8")
    values = [re.sub(r"<[^>]+>", "", v) for v in re.findall(r"<si>(.*?)</si>", shared, re.S)]
    village_no = str(village).split("_")[-1].lstrip("0") or "0"
    out = {}
    for m in re.finditer(r'<row r="(\d+)"[^>]*>(.*?)</row>', sheet, re.S):
        rid, body = int(m.group(1)), m.group(2)
        g = re.search(r'<c r="G%d"[^>]*><v>(\d+)</v></c>' % rid, body)
        h = re.search(r'<c r="H%d"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        j = re.search(r'<c r="J%d"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        if g and h and g.group(1) == village_no:
            out[values[int(h.group(1))]] = values[int(j.group(1))] if j else ""
    return out


def build(village, name, out_dir, workbook=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    status = tracker_status(village, workbook)
    rows = []
    for s, (f, _layer) in sorted(visible.shown(village).items(), key=lambda kv: paths.survey_sort_key(kv[0])):
        g = visible.read_hand(f)
        if g.crs is None or g.crs.to_epsg() != 32644:
            g = g.to_crs(32644)
        g = g[g.geometry.notna() & ~g.geometry.is_empty]
        body = unary_union([x.buffer(0) for x in g.geometry])
        iou = None
        if paths.sheet_path(village, s).exists():
            sheet = sheets.load_sheet(village, s)
            drawing = unary_union([p.buffer(0) for _pr, p in sheet])
            iou = fabric_check.shape_iou(drawing, body)[0]
        for _, r in g.iterrows():
            sub, kide = fmb_finish.survey_keys(s, r.get("plot_no"))
            a = r.geometry.area
            rows.append({"kide": kide, "survey_no": s, "subdiv_no": sub, "area_acre": round(a / ACRE, 4),
                         "area_sqm": round(a, 2), "village_code": village, "village": name,
                         "status": status.get(s, ""), "iou_vs_pdf": iou, "source_file": f.name,
                         "geometry": r.geometry})
    plots = gpd.GeoDataFrame(rows, geometry="geometry", crs=32644).to_crs(4326)
    plots_path = out_dir / ("%s_Manual_Georeferenced.geojson" % name)
    plots.to_file(plots_path, driver="GeoJSON")
    check = fabric_check.check_village(village)
    csv_path = fabric_check.write_csv(village, check, out_dir)
    return {"plots": plots_path, "check": csv_path, "n_plots": len(plots),
            "n_surveys": plots.survey_no.nunique(), "overlaps": len(check["overlaps"]),
            "gaps": len(check["gaps"])}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("village")
    ap.add_argument("name")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    print(build(a.village, a.name, a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
