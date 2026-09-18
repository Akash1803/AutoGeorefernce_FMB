r"""
offpage_check.py
================
Find sketches whose survey boundary runs off the sheet (long strip surveys such
as railway or road land drawn at 1:500 on A0). Such sheets cannot close into
polygons and must be re-fetched from the portal at a coarser scale.

A sheet is 'off-page' when a main boundary segment reaches the page margin band
(the title/footer band or the side frame), i.e. the drawing was clipped.

    python offpage_check.py            -> prints one line per off-page sheet, writes offpage_sheets.csv
"""
import os, sys, csv, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pdfplumber
import fmb_to_geojson as F

D = r"D:\Projects\Tambaram_Chengalpattu"; SK = os.path.join(D, "FMB_Sketches")


def is_offpage(pdf):
    page = pdfplumber.open(pdf).pages[0]
    W, H = float(page.width), float(page.height); k = H / F.A0_H
    segs = [s for s in F.extract_lines(page) if s['layer'] in ('survey_boundary_main', 'subdivision_line')]
    hits = 0
    for s in segs:
        for (x, y) in (s['p0'], s['p1']):          # y already flipped: origin bottom-left
            if y <= 155 * k or y >= H - 155 * k or x <= 60 * k or x >= W - 60 * k:
                hits += 1
    return hits, len(segs)


if __name__ == '__main__':
    rows = []
    for pdf in sorted(glob.glob(os.path.join(SK, "35_*", "*.pdf"))):
        if "_s" in os.path.basename(pdf):
            continue
        hits, n = is_offpage(pdf)
        if hits:
            code = os.path.basename(os.path.dirname(pdf)); survey = os.path.splitext(os.path.basename(pdf))[0]
            rows.append({"village_code": code, "survey_no": survey, "clipped_endpoints": hits, "boundary_segments": n, "pdf": os.path.relpath(pdf, D)})
            print("%s %s: %d clipped endpoints of %d boundary segments" % (code, survey, hits, n))
    with open(os.path.join(SK, "offpage_sheets.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["village_code", "survey_no", "clipped_endpoints", "boundary_segments", "pdf"]); w.writeheader(); w.writerows(rows)
    print("off-page sheets:", len(rows))
