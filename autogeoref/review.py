"""Everything the team reads: the status CSV, the corridor worklist, the tracker and the QGIS group.

The tracker is patched as a zip, never round-tripped through openpyxl, which silently deletes the
x14 extension dataValidations holding the Status (Lists!$A$2:$A$8) and Assigned To (Team!$A$2:$A$31)
dropdowns. The Status column belongs to the team and is never written (user instruction,
2026-09-18: "Don't modify my style, colours and other stuffs, just update the sheet is enough").
"""
import csv
import datetime
import json
import re
import shutil
import zipfile
from pathlib import Path

from . import paths, qgis_bridge

STATUS_COLUMNS = ["survey", "status", "method", "colour", "confidence", "share", "margin",
                  "observable", "boundary_rms_m", "matched_length_m", "support_m", "support_next_m",
                  "n_candidates", "printed_far", "side_ok", "side_bad", "contradictions", "certain",
                  "anchor_partners", "pass", "neighbours", "heading_deg",
                  "shift_m", "sigma_pos_m", "sigma_head_deg", "puvi_reference_m", "puvi_trusted",
                  "puvi_ratio", "anchor_file", "notes", "file", "fp_placed", "fp_final", "run"]
TRACKER_COLUMNS = ["Auto colour", "Auto note", "Auto run"]
COLOURS = {"green": "#19e68c", "amber": "#ffb000", "red": "#ff4136", "anchor": "#000000"}


def write_status(village, rows):
    p = paths.status_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=STATUS_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: paths.survey_sort_key(r["survey"])):
            w.writerow({c: r.get(c, "") for c in STATUS_COLUMNS})
    return p


def read_status(village):
    p = paths.status_path(village)
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def tool_written_fingerprints(village):
    out = set()
    for r in read_status(village):
        for c in ("fp_placed", "fp_final"):
            if r.get(c):
                out.add(r[c])
    return out


def write_worklist(rows_by_village):
    p = paths.PROJECT / "FMB_Vector" / "worklist.csv"
    p.parent.mkdir(parents=True, exist_ok=True)
    cols = ["village_code", "survey", "colour", "confidence", "method", "notes", "file"]
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for village in sorted(rows_by_village):
            for r in sorted(rows_by_village[village], key=lambda r: paths.survey_sort_key(r["survey"])):
                w.writerow(dict({c: r.get(c, "") for c in cols}, village_code=village))
    return p


# --------------------------------------------------------------------------- tracker (XML patch)
def _escape(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _shared_strings(text):
    return ["".join(re.findall(r"<t[^>]*>(.*?)</t>", si, re.S))
            for si in re.findall(r"<si>(.*?)</si>", text, re.S)]


def _next_col(col):
    return chr(ord(col) + 1) if len(col) == 1 and col < "Z" else col + "A"


def _index(col):
    n = 0
    for ch in col:
        n = n * 26 + ord(ch) - 64
    return n


def _widths(sheet, first, widths):
    """Give the tool's own columns a readable width without touching the team's."""
    m = re.search(r'<col min="(\d+)" max="16384"[^>]*/>', sheet)
    if not m:
        return sheet
    style = re.search(r'style="(\d+)"', m.group(0))
    style = ' style="%s"' % style.group(1) if style else ""
    own = "".join('<col min="%d" max="%d" width="%d"%s customWidth="1"/>'
                  % (first + i, first + i, w, style) for i, w in enumerate(widths))
    tail = '<col min="%d" max="16384" width="8.7265625"%s/>' % (first + len(widths), style)
    head = ('<col min="%s" max="%d" width="8.7265625"%s/>' % (m.group(1), first - 1, style)
            if int(m.group(1)) < first else "")
    return sheet.replace(m.group(0), head + own + tail, 1)


def update_tracker(rows, workbook=None):
    """Add the three tool-owned columns, keyed by village code + survey number.

    Returns 'applied', 'locked', 'no change' or 'no workbook'.
    """
    workbook = Path(workbook or paths.TRACKER)
    if not workbook.exists():
        return "no workbook"
    if (workbook.parent / ("~$" + workbook.name)).exists():
        return "locked"
    zin = zipfile.ZipFile(workbook)
    parts = {n: zin.read(n) for n in zin.namelist()}
    zin.close()
    sheet = parts["xl/worksheets/sheet1.xml"].decode("utf-8")
    shared = parts["xl/sharedStrings.xml"].decode("utf-8")
    values = _shared_strings(shared)

    header = re.search(r'<row r="1".*?</row>', sheet, re.S).group(0)
    used = re.findall(r'<c r="([A-Z]+)1"', header)
    last_col = used[-1]
    rowmap = {}
    for m in re.finditer(r'<row r="(\d+)".*?</row>', sheet, re.S):
        rid, body = int(m.group(1)), m.group(0)
        if rid == 1:
            continue
        g = re.search(r'<c r="G%d"[^>]*><v>(\d+)</v></c>' % rid, body)
        h = re.search(r'<c r="H%d"[^>]*t="s"><v>(\d+)</v></c>' % rid, body)
        if g and h:
            rowmap[(g.group(1), values[int(h.group(1))])] = rid

    new_strings = []

    def sref(text):
        if text in values:
            return values.index(text)
        if text in new_strings:
            return len(values) + new_strings.index(text)
        new_strings.append(text)
        return len(values) + len(new_strings) - 1

    cols, cursor = {}, last_col
    for name in TRACKER_COLUMNS:
        cursor = _next_col(cursor)
        cols[name] = cursor
        sheet = sheet.replace(header, header.replace(
            "</row>", '<c r="%s1" s="7" t="s"><v>%d</v></c></row>' % (cursor, sref(name))), 1)
        header = re.search(r'<row r="1".*?</row>', sheet, re.S).group(0)

    changed = 0
    for r in rows:
        village_no = str(r.get("village_code", "")).split("_")[-1].lstrip("0") or "0"
        rid = rowmap.get((village_no, str(r.get("survey", ""))))
        if rid is None:
            continue
        cells = "".join('<c r="%s%d" s="9" t="s"><v>%d</v></c>' % (cols[name], rid, sref(str(val)))
                        for name, val in (("Auto colour", r.get("colour", "")),
                                          ("Auto note", r.get("notes", "")),
                                          ("Auto run", r.get("run", ""))) if str(val))
        if not cells:
            continue
        sheet = re.sub(r'(<row r="%d".*?)</row>' % rid, lambda m: m.group(1) + cells + "</row>",
                       sheet, count=1, flags=re.S)
        changed += 1
    if not changed:
        return "no change"
    if new_strings:
        n0 = int(re.search(r'uniqueCount="(\d+)"', shared).group(1))
        c0 = int(re.search(r'\bcount="(\d+)"', shared).group(1))
        shared = shared.replace("</sst>", "".join("<si><t>%s</t></si>" % _escape(s)
                                                  for s in new_strings) + "</sst>")
        shared = shared.replace('count="%d" uniqueCount="%d"' % (c0, n0),
                                'count="%d" uniqueCount="%d"' % (c0 + changed * 3, n0 + len(new_strings)), 1)
        parts["xl/sharedStrings.xml"] = shared.encode("utf-8")
    last = cols[TRACKER_COLUMNS[-1]]
    sheet = re.sub(r'<dimension ref="A1:[A-Z]+(\d+)"/>',
                   lambda m: '<dimension ref="A1:%s%s"/>' % (last, m.group(1)), sheet, count=1)
    # the filter must reach the new columns, or they read as frozen next to a filtered table
    sheet = re.sub(r'(<autoFilter ref="A1:)[A-Z]+(\d+")',
                   lambda m: m.group(1) + last + m.group(2), sheet, count=1)
    sheet = _widths(sheet, _index(last) - len(TRACKER_COLUMNS) + 1, [10, 40, 20])
    parts["xl/worksheets/sheet1.xml"] = sheet.encode("utf-8")
    parts["xl/workbook.xml"] = re.sub(rb'<calcPr(?![^>]*fullCalcOnLoad)', b'<calcPr fullCalcOnLoad="1"',
                                      parts["xl/workbook.xml"], count=1)
    backup = paths.logs_dir() / "tracker_backups"
    backup.mkdir(parents=True, exist_ok=True)
    shutil.copy2(workbook, backup / (workbook.stem + "_before_autogeoref_%s.xlsx"
                                     % datetime.datetime.now().strftime("%Y%m%d%H%M%S")))
    tmp = workbook.with_suffix(".xlsx.autogeoref-tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for name, data in parts.items():
            zo.writestr(name, data)
    try:
        tmp.replace(workbook)
    except PermissionError:
        tmp.unlink(missing_ok=True)
        return "locked"
    return "applied"


# --------------------------------------------------------------------------- QGIS review group
_GROUP = """
import json, os
from qgis.core import (QgsProject, QgsVectorLayer, QgsFillSymbol, QgsSingleSymbolRenderer,
                       QgsLayerTreeGroup)
proj = QgsProject.instance()
want = %r
if want and os.path.normcase(proj.fileName()) != os.path.normcase(want):
    print(json.dumps({"skipped": "a different project is open", "open": proj.fileName()}))
else:
    root = proj.layerTreeRoot(); name = "Georef review - " + %r
    old = root.findGroup(name)
    if old:
        for n in old.findLayers():
            proj.removeMapLayer(n.layerId())
        root.removeChildNode(old)
    grp = QgsLayerTreeGroup(name); root.insertChildNode(0, grp); added = []
    for path, colour, label in %s:
        lyr = QgsVectorLayer(path + "|layername=parcels", label, "ogr")
        if not lyr.isValid():
            lyr = QgsVectorLayer(path, label, "ogr")
        if not lyr.isValid():
            continue
        lyr.setRenderer(QgsSingleSymbolRenderer(QgsFillSymbol.createSimple(
            {"color": colour + "40", "outline_color": colour, "outline_width": "0.6"})))
        proj.addMapLayer(lyr, False); grp.addLayer(lyr); added.append(label)
    print(json.dumps({"group": name, "layers": added}))
"""


def build_group(village, rows, project_path=None):
    """Load the run's outputs into QGIS. Never fatal: files and CSVs are written regardless."""
    if not qgis_bridge.available():
        return {"skipped": "QGIS not reachable; rerun with --review when QGIS is open"}
    spec = [(str(paths.output_path(village, r["survey"])), COLOURS.get(r.get("colour"), "#888888"),
             "%s (%s)" % (r["survey"], r.get("colour", "")))
            for r in rows if paths.output_path(village, r["survey"]).exists()]
    try:
        return qgis_bridge.json_result(_GROUP % (project_path or "", village, json.dumps(spec)))
    except Exception as exc:
        return {"skipped": str(exc)}
