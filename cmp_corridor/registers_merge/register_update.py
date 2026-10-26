"""Refresh the Tamilnilam + Reginet columns of every plot of Akash's selected surveys, with the SAME matching rules as the
delivered master (Merged_Data/_scripts/03_crosscheck.py + 05_master_excel.py):
 Tamilnilam - village + survey unit (109A) + subdivision exact, status Success; a plot without a number only when the register
              holds the survey undivided.
 Reginet    - glv_of() of 03_crosscheck (exact; parent; children sharing one value; whole survey; ...), copied verbatim.
Writes fix_registers.json {fid-less: kide, vc, pt, attrs incl. None} + registers_report.csv. Changes attributes only."""
import sys, re, json, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd, geopandas as gpd
D = "D:/Projects/CMP_Roads_LPS/Merged_Data"
src = open(f"{D}/_scripts/03_crosscheck.py", encoding="utf-8").read()
ns = {"pd": pd, "re": re}
exec("n = " + src.split("\nn = ", 1)[1].split("\n", 1)[0], ns)
blk = src.split("G = pd.read_excel", 1)[1].split("gb = [glv_of", 1)[0]
exec("G = pd.read_excel" + blk.replace("os.path.join(D, \"04_Reginet_output.xlsx\")", repr(f"{D}/04_Reginet_output.xlsx")), ns)
n, glv_of = ns["n"], ns["glv_of"]
T = pd.read_excel(f"{D}/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
T["vc"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3)
T["unit"] = T.fmb_survey_no.map(n); T["sub"] = T.sub_division.map(n)
reg = {(r.vc, r.unit, r.sub): r for r in T.itertuples()}
units = T.groupby(["vc", "unit"]).sub.apply(set).to_dict()


def split_val(v):
    v = str(v)
    if not v.strip() or v.strip() == "-":
        return None, None
    a, _, b = v.partition("/")
    try:
        return float(a.strip()), b.strip()
    except ValueError:
        return None, None


M = gpd.read_file(os.environ.get("REG_FILE", "D:/Projects/CMP_Roads_LPS/FMB_Georeferenced/CMP_my_Part.geojson"))   # Akash's working layer (2026-10-23 night): every plot of it
M["kide"] = M.kide.astype(str); M["vc"] = M.village_code.astype(str); M["sv"] = M.kide.str.split("/").str[0]
S = M
TN = ["tamilnilam_extent_sqm", "tamilnilam_patta_number", "tamilnilam_land_type", "tamilnilam_government_priority", "tamilnilam_extent_hectares",
      "tamilnilam_extent_ares", "tamilnilam_extent_acre", "tamilnilam_owner_details", "tamilnilam_remarks", "tamilnilam_remarks1"]
RG = ["reginet_guideline_value", "reginet_unit", "reginet_guideline_value_metric", "reginet_unit_metric", "reginet_land_classification"]
same = lambda a, b: (None if a is None or (isinstance(a, float) and a != a) or str(a) == "" else str(a)) == (None if b is None or (isinstance(b, float) and b != b) or str(b) == "" else str(b)) or \
    (_num(a) is not None and _num(a) == _num(b))


def _num(x):
    try:
        return round(float(x), 4)
    except (TypeError, ValueError):
        return None


out, rep = [], []
for r in S.itertuples():
    unit, sub = n(r.sv), n(r.subdiv_no)
    subs = units.get((r.vc, unit)); t = None
    if subs is None:
        tm = "survey not in the register"
    elif sub:
        t = reg.get((r.vc, unit, sub)); tm = "matched" if t is not None else "subdivision not in the register"
    elif subs <= {""}:
        t = reg.get((r.vc, unit, "")); tm = "matched (undivided survey)"
    else:
        tm = f"plot has no subdivision number; register has {len(subs - {''})} subdivisions"
    if t is not None and t.tamilnilam_status != "Success":
        t, tm = None, "Tamilnilam returned no data"
    new = {}
    if t is not None:
        sqm = float(t.extent_hectares or 0) * 10000 + float(t.extent_ares or 0) * 100
        new.update(dict(zip(TN, [round(sqm, 2), t.patta_number, t.land_type_english, t.government_priority, t.extent_hectares, t.extent_ares,
                                 round(sqm / 4046.856, 6), t.owner_details, t.remarks, t.remarks1])))
    else:
        new.update({c: None for c in TN})
    g, basis = glv_of(r.vc, unit, sub, new["tamilnilam_land_type"] or "")
    vb, ub = split_val(g["Guideline_Value_British"]) if g is not None else (None, None)
    vm, um = split_val(g["Guideline_Value_Metric"]) if g is not None else (None, None)
    if g is not None and vm is not None:          # as 05_master_excel: classification only with a value; "0" (Govt.) is a value
        new.update(dict(zip(RG, [vb, ub, vm, um, g["Land_Classification"]])))
    else:                                    # no rate on Reginet ("-", "Refer Street Value" ...): all Reginet fields left empty (Akash 2026-10-23)
        new.update({c: None for c in RG}); basis = "no guideline value on Reginet"
    diff = {k: v for k, v in new.items() if not same(v, getattr(r, k))}
    rep.append({"village_code": r.vc, "kide": r.kide, "tamilnilam": tm, "reginet": basis or "matched", "reginet_subdivision": g["Subdivision Number"] if g is not None else "",
                "changed_fields": ", ".join(diff), "old_extent_sqm": r.tamilnilam_extent_sqm, "new_extent_sqm": new["tamilnilam_extent_sqm"],
                "old_glv": r.reginet_guideline_value, "new_glv": new["reginet_guideline_value"], "plot_area_sqm": r.area_sqm})
    if diff:
        rp = r.geometry.representative_point()
        out.append({"vc": r.vc, "old_kide": r.kide, "old_pt": [rp.x, rp.y], "attrs": diff})
R = pd.DataFrame(rep); R.to_csv(os.environ.get("REG_OUT", "fix_registers.json").replace(".json", "_report.csv"), index=False, encoding="utf-8-sig")
json.dump({"change": out}, open(os.environ.get("REG_OUT", "fix_registers.json"), "w"), default=str)
print("plots", len(S), "| changed", len(out))
print(R.tamilnilam.value_counts().to_dict()); print(R.reginet.str.split(" ").str[:2].str.join(" ").value_counts().to_dict())
