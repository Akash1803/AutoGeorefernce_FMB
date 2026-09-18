r"""
refresh_sketch_summary.py
=========================
Rebuild FMB_Sketches\fmb_sketch_summary.csv and fmb_missing_sketches.csv from
what is actually on disk:

* base survey numbers  : _logs\fmb_batch_input.csv (one row per village + numeric survey number)
* sketches             : FMB_Sketches\<dd_tt_vvv>\<unit>.pdf  (unit = plain number or letter form)
* portal error reasons : _logs\download_log.csv

A base number is "covered" by its plain-number sketch, or by its letter-form
sketches (64A, 64B ...) when the plain number does not exist on the portal.
"""
import os, re, csv, glob, collections
D = r"D:\Projects\Tambaram_Chengalpattu"; SK = os.path.join(D, "FMB_Sketches"); LOGS = os.path.join(D, "_logs")
base_rows = list(csv.DictReader(open(os.path.join(LOGS, "fmb_batch_input.csv"), encoding="utf-8")))
log = collections.defaultdict(list)
for e in csv.DictReader(open(os.path.join(LOGS, "download_log.csv"), encoding="utf-8")):
    log[e["giscode"]].append(e)
def last_error(gis):
    es = [e for e in log.get(gis, []) if e["status"] != "saved"]
    return es[-1]["api_error"] if es else ""
on_disk = collections.defaultdict(set)
for p in glob.glob(os.path.join(SK, "35_*", "*.pdf")):
    if os.path.getsize(p) > 1000:
        on_disk[os.path.basename(os.path.dirname(p))].add(os.path.splitext(os.path.basename(p))[0])
F = ["village", "village_code", "survey_unit", "base_survey", "unit_kind", "sketch", "pdf", "api_error", "vector_pieces"]
summary, missing = [], []
claimed = collections.defaultdict(set)
for r in base_rows:
    t, v = r["taluk"].zfill(2), r["village"].zfill(3); code = "35_%s_%s" % (t, v); num = r["survey_no"]
    if num in on_disk[code]:
        summary.append(dict(village=r["village_name"], village_code=code, survey_unit=num, base_survey=num, unit_kind="plain number",
                            sketch="saved", pdf="%s\\%s.pdf" % (code, num), api_error="", vector_pieces=r["vector_pieces"]))
        claimed[code].add(num); continue
    letters = sorted(u for u in on_disk[code] if re.fullmatch(num + r"[A-Za-z]+", u))
    for u in letters:
        summary.append(dict(village=r["village_name"], village_code=code, survey_unit=u, base_survey=num, unit_kind="letter form",
                            sketch="saved", pdf="%s\\%s.pdf" % (code, u), api_error="", vector_pieces=r["vector_pieces"]))
        claimed[code].add(u)
    if not letters:
        pieces = [re.sub(r"\s+", "", p) for p in r["vector_pieces"].split(";")]
        tried = sorted({p for p in pieces if re.search(r"[A-Za-z]", p)} | {num + s for s in "ABCDEF" if "S35%s%s%s%s" % (t, v, num, s) in log})
        row = dict(village=r["village_name"], village_code=code, survey_unit=num, base_survey=num, unit_kind="plain number", sketch="missing", pdf="",
                   api_error=(last_error("S35%s%s%s" % (t, v, num)) or "not available") + " | also tried: " + (", ".join(tried) or "-"), vector_pieces=r["vector_pieces"])
        summary.append(row); missing.append(row)
# sketches on disk that no base number claimed (e.g. added by hand)
for code, units in on_disk.items():
    for u in sorted(units - claimed[code]):
        vname = next((r["village_name"] for r in base_rows if "35_%s_%s" % (r["taluk"].zfill(2), r["village"].zfill(3)) == code), code)
        summary.append(dict(village=vname, village_code=code, survey_unit=u, base_survey=re.match(r"^\d+", u).group(0) if re.match(r"^\d+", u) else u,
                            unit_kind="letter form" if re.search(r"[A-Za-z]", u) else "plain number", sketch="saved", pdf="%s\\%s.pdf" % (code, u), api_error="", vector_pieces="(added manually)"))
for name, rows in (("fmb_sketch_summary.csv", summary), ("fmb_missing_sketches.csv", missing)):
    try:
        with open(os.path.join(SK, name), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=F); w.writeheader(); w.writerows(rows)
    except PermissionError:
        print("LOCKED (open in Excel?) - not rewritten:", name)
saved = [s for s in summary if s["sketch"] == "saved"]
print("survey units:", len(summary), "| sketches:", len(saved), "(plain %d, letter forms %d)" % (sum(1 for s in saved if s["unit_kind"] == "plain number"), sum(1 for s in saved if s["unit_kind"] == "letter form")),
      "| base numbers with no sketch:", len({(m["village_code"], m["base_survey"]) for m in missing}))
print("still missing:", collections.Counter(m["village"] for m in missing).most_common())
