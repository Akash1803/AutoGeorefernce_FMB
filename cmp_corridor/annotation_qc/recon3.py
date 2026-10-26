"""Three-way reconciliation per LAYER plot: machine read (ann_check, label-anchored warp) + two independent tile readers.
Tile reads are joined to the machine read's PDF plots by centroid (qc_centroids.csv vs <pfx>_pdfplots.csv, < 10 pt), never by
the tiles' own rigid layer attribution.
Status: OK (all reads = label) | FIX (machine + >= 1 reader agree on another number, or both readers agree) | CHECK (any
disagreement, single-read wrongs, unmatched plots, missing numbers nobody confirms) | NOT READ.
Usage: python recon3.py <pfx> <tiledir> [<tiledir>...] -> qc/<pfx>_decisions.csv"""
import sys, glob, csv, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
pfx, tiledirs = sys.argv[1], sys.argv[2:]
norm = lambda s: re.sub(r"[^0-9A-Z]", "", str(s).upper()) if str(s).strip().lower() not in ("", "-", "nan", "none", "?", "null") else ""
M = pd.read_csv(f"qc/{pfx}_read1.csv", dtype=str).fillna(""); P = pd.read_csv(f"qc/{pfx}_pdfplots.csv", dtype=str).fillna("")
M["pdf_plot"] = M.pdf_plot.str.replace(r"\.0$", "", regex=True)
_T = pd.read_excel("D:/Projects/CMP_Roads_LPS/Merged_Data/02_Tamilnilam_output.xlsx", dtype=str).fillna("")
_T["vc"] = _T.district_id.str.zfill(2) + "_" + _T.taluk_id.str.zfill(2) + "_" + _T.village_id.str.zfill(3)
REG = _T.groupby([_T.vc, _T.fmb_survey_no.map(norm)]).sub_division.apply(lambda x: {norm(v) for v in x} - {""}).to_dict(); P["pdf_plot"] = P.pdf_plot.str.replace(r"\.0$", "", regex=True)
M = M[(M.kide != "") & ~M.verdict.isin(["LABEL NOT IN REGISTER", "PRINTED NUMBER NOT IN LAYER"]) & ~M.verdict.str.startswith("DUPLICATE")].copy(); M["key"] = M.village + "|" + M.kide
# qc_id -> (village, survey, pdf_plot)
qmap = {}
for td in tiledirs:
    C = pd.read_csv(f"{td}/qc_centroids.csv", dtype=str)
    for (v, s), g in C.groupby(["village", "survey"]):
        p = P[(P.village == v) & (P.survey == s)]
        if not len(p): continue
        px = p[["x_pt", "y_pt"]].astype(float).values
        for r in g.itertuples():
            d = np.hypot(*(px - np.array([float(r.x_pt), float(r.y_pt)])).T)
            if d.min() < 10: qmap[(td, v, r.qc_id)] = (v, s, p.pdf_plot.iloc[int(d.argmin())])
reads = []
import os
for td in tiledirs:
    Qp = pd.read_csv(f"{td}/plots.csv", dtype=str).fillna("")
    for f in glob.glob(f"{td}/reads/*.csv"):
        fn = f.replace("\\", "/").split("/")[-1]; who = fn[0]; grp = fn.split("_")[1].split(".")[0]
        gfile = f"{td}/{grp}.txt"; tiles = open(gfile).read().split() if os.path.exists(gfile) else []
        vill_of = {}                                   # survey -> villages whose tiles are in this group
        for t_ in tiles:
            base = t_.rsplit("_t", 1)[0]; vc_, sv_ = base[:9], base[10:]; vill_of.setdefault(sv_, set()).add(vc_)
        for r in list(csv.reader(open(f, encoding="utf-8")))[1:]:
            if len(r) < 3: continue
            qid = r[0].strip(); sv_ = qid.split(".")[0]; cands = sorted(vill_of.get(sv_, set()))
            if len(cands) > 1:                         # same survey number in two villages inside this group
                note = r[4] if len(r) > 4 else ""
                byn = [c for c in cands if c in note]
                if len(byn) == 1: cands = byn
                else:
                    seen = norm(r[1]); byl = [c for c in cands if any(norm(x) == seen for x in Qp[(Qp.village == c) & (Qp.qc_id == qid)].layer_label)]
                    cands = byl if len(byl) == 1 else []
            if len(cands) != 1: continue
            vc_ = cands[0]
            if (td, vc_, qid) not in qmap: continue
            v, s, k = qmap[(td, vc_, qid)]; verdict = r[2].strip().lower(); printed = r[3].strip() if len(r) > 3 else ""
            read = None if verdict == "unsure" else (norm(r[1]) if verdict == "ok" else "" if verdict == "none" else norm(printed))
            reads.append({"village": v, "survey": s, "pdf_plot": k, "who": who, "read": read})
T = pd.DataFrame(reads); print("tile reads joined:", len(T), "| qmap:", len(qmap), "| by dir:", {td: sum(1 for k in qmap if k[0] == td) for td in tiledirs})
dec = []
for r in M.itertuples():
    lab = norm(r.layer_label)
    mnum = sorted({norm(t.rstrip("*")) for t in r.pdf_numbers.split("|") if t})
    if r.verdict in ("UNMATCHED", "NO PDF", "NO PDF PLOTS", "LABEL NOT IN REGISTER", "PRINTED NUMBER NOT IN LAYER") or r.verdict.startswith("DUPLICATE"): m = None
    elif r.verdict == "ok (no number)": m = ""
    elif lab and lab in mnum: m = lab
    elif len(mnum) == 1: m = mnum[0]
    else: m = None
    t = T[(T.village == r.village) & (T.survey == r.survey) & (T.pdf_plot == r.pdf_plot)] if r.pdf_plot else T.iloc[0:0]
    rd = {}
    for w in ("A", "B"):
        x = t[t.who == w]
        if len(x): rd[w] = x.read.iloc[0]
    said = {k: v for k, v in {"M": m, **rd}.items() if v is not None}
    regset = REG.get((r.village, norm(r.survey)), set())
    readers = [v for k, v in said.items() if k != "M"]
    new = ""
    both_readers = len(readers) == 2 and readers[0] == readers[1]
    if not lab and both_readers and readers[0] == "": status = "OK"                                        # unnumbered piece, both readers: nothing printed
    elif r.verdict in ("UNMATCHED", "NO PDF", "NO PDF PLOTS"):
        if both_readers and readers[0] == lab and lab: status = "OK"
        elif both_readers and readers[0] not in ("", lab) and (readers[0] in regset or not regset): status = "FIX (readers)"; new = readers[0]
        else: status = "CHECK (no PDF match)"
    elif not said: status = "NOT READ"
    elif all(v == lab for v in said.values()): status = "OK"
    elif m is not None and m != lab and m != "" and m in readers: status = "FIX"; new = m                      # machine + a reader agree
    elif m is None and len(readers) == 2 and readers[0] == readers[1] and readers[0] not in ("", lab) and (readers[0] in regset or not regset): status = "FIX (readers)"; new = readers[0]
    elif len(readers) == 2 and readers[0] == readers[1] == "" and lab: status = "CHECK (reads say no number printed)"
    else: status = "CHECK"
    dec.append({"key": r.key, "village": r.village, "survey": r.survey, "kide": r.kide, "layer_label": r.layer_label, "machine": "?" if m is None else (m or "(none)"),
                "machine_verdict": r.verdict, "A": "-" if "A" not in rd else ("unsure" if rd["A"] is None else (rd["A"] or "(none)")),
                "B": "-" if "B" not in rd else ("unsure" if rd["B"] is None else (rd["B"] or "(none)")), "status": status, "new_label": new, "iou": r.iou, "pdf_plot": r.pdf_plot})
D = pd.DataFrame(dec); D.to_csv(f"qc/{pfx}_decisions.csv", index=False, encoding="utf-8-sig")
print(D.status.value_counts().to_dict())
pd.set_option("display.width", 230); pd.set_option("display.max_rows", 500)
print(D[D.status.str.startswith("FIX")][["village", "kide", "layer_label", "machine", "A", "B", "status", "new_label", "iou"]].to_string(index=False))
