r"""
batch_fmb_to_vector.py
======================
Run fmb_to_geojson.run() over every downloaded FMB sketch of the Tambaram -
Chengalpattu rail corridor and write a summary.

    python batch_fmb_to_vector.py [--only 35_05_134] [--limit N]

Inputs : D:\Projects\Tambaram_Chengalpattu\FMB_Sketches\<dd_tt_vvv>\<survey>.pdf  (+ sheet_scales.csv)
Outputs: D:\Projects\Tambaram_Chengalpattu\FMB_Vector\<dd_tt_vvv>\<survey>_parcels.geojson
         (+ _parcels.csv, _lines.geojson, _preview.png per sketch)
         D:\Projects\Tambaram_Chengalpattu\FMB_Vector\fmb_vector_summary.csv

Scale: sheets are fetched at 1:500; long strip surveys that ran off the page were
re-fetched at a coarser scale (recorded per sheet in FMB_Sketches\sheet_scales.csv,
file names stay plain). The sheet's own edge dimensions verify the scale.
"""
import os, sys, csv, glob, json, time, argparse, contextlib, io, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fmb_to_geojson as F
import scale_from_dimensions as SD

D = r"D:\Projects\Tambaram_Chengalpattu"
SK = os.path.join(D, "FMB_Sketches"); OUT = os.path.join(D, "FMB_Vector")
GLYPHS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "glyphlib.json")
LIB = [(c, tuple(b), x) for c, b, x in json.load(open(GLYPHS))]
TALUK = {"04": "Chengalpattu", "05": "Tambaram", "15": "Vandalur"}


def survey_areas():
    """(village_code, survey_no) -> total plot area (m2) from the village vector layers."""
    import sqlite3, re
    con = sqlite3.connect(os.path.join(D, "Railway_Buffer_Vector_Plots.gpkg"))
    out = collections.defaultdict(float); names = {}
    for code, s, a, v in con.execute("select village_code, survey_no, plot_area_sqm, village from vector_in_buffer_30m"):
        m = re.match(r"^\d+", s.strip())
        if m:
            out[(code, m.group(0))] += a or 0; names[code] = v
    return out, names


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--only"); ap.add_argument("--limit", type=int)
    ap.add_argument("--sheets", help="comma list of village_code/survey to (re)run, e.g. 35_04_074/13,35_15_002/224A")
    ap.add_argument("--merge", action="store_true", help="with --sheets: update those rows in the existing summary instead of rewriting it")
    a = ap.parse_args()
    areas, names = survey_areas()
    import re
    pick = {}
    for pth in sorted(glob.glob(os.path.join(SK, "35_*", "*.pdf"))):
        stem = os.path.splitext(os.path.basename(pth))[0]
        m = re.match(r"^(.+?)(?:_s(\d+))?$", stem); base, sc = m.group(1), int(m.group(2) or 500)
        key = (os.path.basename(os.path.dirname(pth)), base)
        if key not in pick or sc > pick[key][1]:
            pick[key] = (pth, sc)
    pdfs = [v[0] for v in sorted(pick.values())]
    nominal = {v[0]: v[1] for v in pick.values()}
    # sheets re-fetched at another scale keep their plain file name; the scale lives in sheet_scales.csv
    sc_file = os.path.join(SK, "sheet_scales.csv")
    if os.path.exists(sc_file):
        for r in csv.DictReader(open(sc_file, encoding="utf-8")):
            pth = os.path.join(SK, r["village_code"], r["survey_no"] + ".pdf")
            if pth in nominal:
                nominal[pth] = int(r["scale"])
    if a.only: pdfs = [p for p in pdfs if os.path.basename(os.path.dirname(p)) == a.only]
    if a.limit: pdfs = pdfs[:a.limit]
    if a.sheets:
        want = set(a.sheets.split(","))
        pdfs = [p for p in pdfs if "%s/%s" % (os.path.basename(os.path.dirname(p)), re.sub(r"_s\d+$", "", os.path.splitext(os.path.basename(p))[0])) in want]
    os.makedirs(OUT, exist_ok=True)
    rows = []; t0 = time.time()
    for i, pdf in enumerate(pdfs, 1):
        code = os.path.basename(os.path.dirname(pdf)); survey = re.sub(r"_s\d+$", "", os.path.splitext(os.path.basename(pdf))[0])
        nom = nominal.get(pdf, 500)
        vname = names.get(code, code); taluk = TALUK.get(code.split("_")[1], "")
        out_dir = os.path.join(OUT, code); vec_area = areas.get((code, re.sub(r"[^0-9].*$", "", survey)))
        row = {"village": vname, "village_code": code, "taluk": taluk, "survey_no": survey, "pdf": os.path.relpath(pdf, D),
               "scale": 500, "scale_note": "", "scale_confidence": "", "polygons": 0, "numbered": 0, "review": 0, "unplaced_labels": 0,
               "polygon_total_sqm": None, "vector_area_sqm": round(vec_area, 1) if vec_area else None, "area_ratio": None,
               "max_node_shift_m": None, "components": None, "text_class": "", "status": "", "error": "",
               "geojson": os.path.relpath(os.path.join(out_dir, os.path.splitext(os.path.basename(pdf))[0] + "_parcels.geojson"), D)}
        try:
            # true scale from the sheet's own edge dimensions (annotated metres vs drawn length)
            # The portal honours the requested scale, so the sheet IS at `nom`; the
            # dimension check only verifies it (annotated metres vs drawn length).
            est = SD.estimate(pdf, LIB, nom); scale = nom; r = est.get("median_ratio")
            if est["confidence"] == "low":
                verdict = "unverified (%d dimensions read)" % est["n"]
            elif r is not None and abs(r - 1.0) <= 0.06:
                verdict = "verified (n=%d, ratio=%.3f)" % (est["n"], r)
            else:
                verdict = "MISMATCH - check manually (n=%d, ratio=%s)" % (est["n"], r)
            row["scale_note"] = ("sheet fetched at 1:%d (strip survey off-page at 1:500); " % nom if nom != 500 else "") + verdict
            with contextlib.redirect_stdout(io.StringIO()):
                feats = F.run(pdf, scale, out_dir, kind="sketch", glyphlib=GLYPHS, area=vec_area,
                              village="%s [%s]" % (vname, code), taluk=taluk, survey=survey)
            st = F.run.last_stats
            total = st["total_sqm"]
            row["scale_confidence"] = est["confidence"]
            row.update({"scale": scale, "polygons": st["polygons"], "numbered": st["numbered"], "review": st["review"],
                        "unplaced_labels": st["unplaced_labels"], "polygon_total_sqm": total,
                        "area_ratio": round(total / vec_area, 3) if vec_area else None,
                        "max_node_shift_m": st["max_node_shift_m"], "components": st["components"], "text_class": st["text_class"],
                        "status": "ok" if st["polygons"] else "no-polygons"})
        except Exception as e:
            row.update({"status": "error", "error": repr(e)[:200]})
        rows.append(row)
        print("[%3d/%d] %s %s: %s poly=%s numbered=%s ratio=%s %s" % (i, len(pdfs), code, survey, row["status"], row["polygons"], row["numbered"], row["area_ratio"], row["scale_note"]), flush=True)
    summ = os.path.join(OUT, "fmb_vector_summary.csv")
    if a.merge and os.path.exists(summ):
        old = list(csv.DictReader(open(summ, encoding="utf-8"))); new = {(r["village_code"], r["survey_no"]): r for r in rows}
        rows = [new.pop((r["village_code"], r["survey_no"]), r) for r in old] + list(new.values())
    with open(summ, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    ok = [r for r in rows if r["status"] == "ok"]
    for r in ok:
        r["polygons"] = int(r["polygons"]); r["numbered"] = int(r["numbered"])   # merged rows come back from CSV as text
    print("\nDONE %d sheets in %.0f s | ok %d | no-polygons %d | errors %d | polygons %d | numbered %d (%.0f%%)" % (
        len(rows), time.time() - t0, len(ok), sum(1 for r in rows if r["status"] == "no-polygons"), sum(1 for r in rows if r["status"] == "error"),
        sum(r["polygons"] for r in ok), sum(r["numbered"] for r in ok), 100 * sum(r["numbered"] for r in ok) / max(1, sum(r["polygons"] for r in ok))))


if __name__ == "__main__":
    main()
