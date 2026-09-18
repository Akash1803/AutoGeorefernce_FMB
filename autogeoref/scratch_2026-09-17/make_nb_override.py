"""Convert the two-reader transcription (merged.json from the read-fmb-neighbours workflow) into nb_override_<village>.json
for autogeoref_demo.sheet_neighbours, and print a comparison with the glyph reader."""
import sys, os, json, re, collections
S = os.path.dirname(os.path.abspath(__file__)); VILLAGE = sys.argv[1]; src = sys.argv[2]
merged = json.load(open(src, encoding="utf-8"))
def norm(num):
    n = str(num).upper().replace(" ", "").replace("S.NO.", "").replace("S.NO", "").replace("SNO", "").split("/")[0]
    return n if re.match(r"^\d+[A-Z]?$", n) else None
out = {}; report = []
for m in merged:
    sv = m["survey"]; rows = []
    for n in m["neighbours"]:
        k = norm(n["number"])
        if not k or k == sv.upper(): continue
        side = collections.Counter(s.upper() for s in n["sides"] if s).most_common(1); side = side[0][0] if side else ""
        conf = "high" if "high" in n["confidences"] and n["readers"] >= 2 else ("medium" if n["readers"] >= 2 or "high" in n["confidences"] else "low")
        rows.append({"number": k, "side": side, "readers": n["readers"], "confidence": conf})
    out[sv] = rows
    report.append((sv, [(r["number"], r["side"], r["readers"]) for r in rows], m.get("railway_side", []), m.get("other_labels", [])))
json.dump(out, open(os.path.join(S, "nb_override_%s.json" % VILLAGE), "w", encoding="utf-8"), indent=1)
glyph = json.load(open(os.path.join(S, "nb_%s.json" % VILLAGE), encoding="utf-8")) if os.path.exists(os.path.join(S, "nb_%s.json" % VILLAGE)) else {}
print("%-5s %-45s %-18s %s" % ("sheet", "readers (number, side, n_readers)", "railway side", "glyph reader"))
for sv, rows, rs, other in report: print("%-5s %-45s %-18s %s | other: %s" % (sv, rows, rs, glyph.get(sv), other[:4]))
print("written", os.path.join(S, "nb_override_%s.json" % VILLAGE))
