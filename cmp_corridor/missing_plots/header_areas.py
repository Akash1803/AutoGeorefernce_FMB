"""Printed survey area (title block) of every FMB sheet under FMB_Sketches -> pdf_header_areas.csv (village_code, survey, pdf_area_m2)."""
import sys, os, glob, csv, time
sys.path.insert(0, "D:/code/FMB_to_GeoJSON")
from autogeoref.printed_numbers import read_header_area
rows = []; done = set()
OUT = "pdf_header_areas.csv"
if os.path.exists(OUT):
    for r in csv.DictReader(open(OUT, encoding="utf-8")): rows.append(r); done.add((r["village_code"], r["survey"]))
files = sorted(glob.glob("D:/Projects/CMP_Roads_LPS/FMB_Sketches/35_*/*.pdf"))
for i, f in enumerate(files):
    vc = os.path.basename(os.path.dirname(f)); sv = os.path.splitext(os.path.basename(f))[0]
    if (vc, sv) in done: continue
    try: a = read_header_area(f)
    except Exception as e: a = None
    rows.append({"village_code": vc, "survey": sv, "pdf_area_m2": "" if a is None else a}); done.add((vc, sv))
    if i % 25 == 0:
        with open(OUT, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["village_code", "survey", "pdf_area_m2"]); w.writeheader(); w.writerows(rows)
with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["village_code", "survey", "pdf_area_m2"]); w.writeheader(); w.writerows(rows)
print("sheets:", len(files), "read:", sum(1 for r in rows if r["pdf_area_m2"] not in ("", None)))
