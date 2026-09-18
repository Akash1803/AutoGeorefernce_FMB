r"""Re-fetch the sheets listed in FMB_Sketches\offpage_sheets.csv from the Collabland
portal at 1:2000 (A0) so that long strip surveys fit the page. Saved beside the
original as <survey>_s2000.pdf. Same endpoint and session as fetch_fmb_gui.py."""
import os, csv, base64, time, requests
D = r"D:\Projects\Tambaram_Chengalpattu\FMB_Sketches"
URL = "https://collabland-tn.gov.in/APIServices/rest/Collabland/FMBMapServicePDF"
s = requests.Session(); s.headers.update({"User-Agent": "Mozilla/5.0", "Accept": "*/*", "Origin": "https://collabland-tn.gov.in",
    "Referer": "https://collabland-tn.gov.in/APIServices/FMBMapService.jsp", "X-Requested-With": "XMLHttpRequest"})
s.cookies.update({"JSESSIONID": "REDACTED_SESSION_COOKIE"})
rows = list(csv.DictReader(open(os.path.join(D, "offpage_sheets.csv"), encoding="utf-8")))
for r in rows:
    _, t, v = r["village_code"].split("_"); gis = "S35%s%s%s" % (t, v, r["survey_no"])
    out = os.path.join(D, r["village_code"], r["survey_no"] + "_s2000.pdf")
    if os.path.exists(out) and os.path.getsize(out) > 1000:
        continue
    try:
        resp = s.post(URL, data={"state": "33", "giscode": gis, "plotno": "", "scale": "2000", "width": "841", "height": "1189", "localLang": "false"}, timeout=90)
        j = resp.json(); b64 = j.get("success") or j.get("data") or j.get("pdf") or j.get("result")
        if b64:
            open(out, "wb").write(base64.b64decode(b64)); print(gis, "-> saved", os.path.basename(out), os.path.getsize(out), "bytes")
        else:
            print(gis, "-> portal error:", str(j.get("error", j))[:80])
    except Exception as e:
        print(gis, "-> failed:", repr(e)[:80])
    time.sleep(1.0)
