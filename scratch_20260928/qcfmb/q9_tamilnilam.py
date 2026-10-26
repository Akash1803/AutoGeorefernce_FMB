"""Strict Tamilnilam match for a table of plots. A plot receives the register record only when
  - its number is verified (both reads, or second read + register, or checked on the PDF), or it is a whole survey;
  - district code, taluk code, village code, survey number and subdivision number are all EQUAL, letter for letter,
    to one record of the A-register that Tamilnilam returned with "Success" under that same survey number;
  - the register extent agrees with the plot (within x2, or within 50 m2 - register extents are rounded to half an are).
No nearest match, no base-number fill, no fuzzy compare. Everything else stays empty and says why.
usage: import q9_tamilnilam as q9; F = q9.match(F)   (F needs village_code, village, survey_no, subdiv_no, area_sqm, subdiv_status)"""
import os, re, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
R = r"D:\Projects\Tambaram_Chengalpattu"
FIELDS = ["patta_number", "land_type_english", "government_priority", "extent_hectares", "extent_ares", "owner_details", "owner_number", "poramboke", "cultivable", "remarks", "remarks1"]
NAMES = {"patta_number": "tamilnilam_patta_number", "land_type_english": "tamilnilam_land_type", "government_priority": "tamilnilam_government_priority", "extent_hectares": "tamilnilam_extent_hectares", "extent_ares": "tamilnilam_extent_ares",
         "owner_details": "tamilnilam_owner_details", "owner_number": "tamilnilam_owner_count", "poramboke": "tamilnilam_poramboke_code", "cultivable": "tamilnilam_cultivable", "remarks": "tamilnilam_remarks", "remarks1": "tamilnilam_remarks1"}
UNVERIFIED = ("to be checked", "not checked", "sliver")
norm = lambda x: re.sub(r"\s+", "", "" if x is None else str(x)).upper()


def child(n, t):
    """n is a subdivision of t: t followed by more, the first added character of the other kind (10 -> 10A -> 10A1; 105 is not under 10)"""
    return len(n) > len(t) and n.startswith(t) and n[len(t)].isdigit() != t[-1].isdigit()


def register():
    T = pd.read_excel(os.path.join(R, "Merged_Data", "02_Tamilnilam_output.xlsx"), dtype=str).fillna("")
    T["village_code"] = T.district_id.str.zfill(2) + "_" + T.taluk_id.str.zfill(2) + "_" + T.village_id.str.zfill(3); T["sqm"] = (T.extent_hectares.replace("", "0").astype(float) * 10000 + T.extent_ares.replace("", "0").astype(float) * 100).round(1)
    T["s"] = T.fmb_survey_no.map(norm); T["n"] = T.sub_division.map(norm); T["own"] = (T.fmb_survey_no == T.clip_survey_no) & (T.tamilnilam_status == "Success")
    assert not T.duplicated(["village_code", "s", "n"]).any()
    return T


def match(F, T=None):
    T = register() if T is None else T; F = F.copy(); ok = T[T.own].set_index(["village_code", "s", "n"]); anyrow = T.set_index(["village_code", "s", "n"])
    subs = {k: sorted(g.n) for k, g in T[T.own].groupby(["village_code", "s"])}; other = {k: g.clip_survey_no.iloc[0] for k, g in T[T.fmb_survey_no != T.clip_survey_no].groupby(["village_code", "s"])}
    F["k_s"] = F.survey_no.map(norm); F["k_n"] = [norm(x) if norm(x) else "-" for x in F.subdiv_no]; F["k_v"] = ~F.subdiv_status.astype(str).str.startswith(UNVERIFIED)
    cnt = F[F.k_v].groupby(["village_code", "k_s", "k_n"]).size().to_dict(); tot = F[F.k_v].groupby(["village_code", "k_s", "k_n"]).area_sqm.sum().to_dict()
    res = []; out = {v: [] for v in NAMES.values()}; ext = []
    for r in F.itertuples():
        k = (r.village_code, r.k_s, r.k_n); rec = None; e = ""
        if not r.k_v: m = "not appended: the plot number is not verified yet"
        elif r.k_n == "-" and not str(r.subdiv_status).startswith(("whole survey", "no FMB sheet")): m = "not appended: the plot has no subdivision number on the sheet"
        elif (r.village_code, r.k_s) not in subs: m = "not appended: the A-register has this survey under the number %s" % other[(r.village_code, r.k_s)] if (r.village_code, r.k_s) in other else "not appended: survey not found in the A-register"
        elif k not in ok.index:
            if k in anyrow.index: m = "not appended: Tamilnilam returned no data for this subdivision (%s)" % anyrow.loc[k].tamilnilam_status
            else:
                S = subs[(r.village_code, r.k_s)]; kids = [t for t in S if t != "-" and r.k_n != "-" and child(t, r.k_n)]; par = [t for t in S if t != "-" and child(r.k_n, t)]
                m = "not appended: " + ("the A-register has the survey in subdivisions (%d)" % len(S) if r.k_n == "-" else "the A-register has it subdivided further (%s)" % ", ".join(kids[:6]) if kids else "the A-register has only the parent plot %s" % max(par, key=len) if par else "subdivision not in the A-register")
        else:
            t = ok.loc[k]; a = float(tot[k]); e = t.sqm; n = cnt[k]
            if t.sqm > 0 and not (0.5 <= a / t.sqm <= 2.0 or abs(a - t.sqm) <= 50): m = "not appended: the register extent differs from the plot (register %.0f m2, plot %.0f m2)" % (t.sqm, a)
            else: rec = t; m = "matched" + (" (subdivision drawn in %d pieces)" % n if n > 1 else "")
        res.append(m); ext.append(e)
        for c, v in NAMES.items(): out[v].append("" if rec is None else rec[c])
    F["tamilnilam_match"] = res; F["tamilnilam_extent_sqm"] = ext
    for v in NAMES.values(): F[v] = out[v]
    F["uid_tn_avail"] = [1 if m.startswith("matched") else 0 for m in res]
    return F.drop(columns=["k_s", "k_n", "k_v"])


def summary(F, name):
    g = F.tamilnilam_match.str.replace(r" \(.*", "", regex=True); print(name, "| plots", len(F), "| matched", int((F.uid_tn_avail == 1).sum())); print(g.value_counts().to_string())
    return pd.crosstab(F.village, g)
