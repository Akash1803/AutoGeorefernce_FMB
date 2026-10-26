"""Dual QC, step 7: what the two readers of the PDF crops saw, set against the machine reads.
A number is taken only when two sources agree: both readers, or one reader and the second (machine) read.
-> out/visual_read.csv (row, final, how) for layer rows; out/visual_read2.csv for the pieces of split rows"""
import os, re, glob, json, warnings; warnings.filterwarnings("ignore")
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
def load(reader):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "out", "visual", f"reader{reader}_batch_*.json"))):
        try: J = json.load(open(f, encoding="utf-8"))
        except Exception as e: print("unreadable", os.path.basename(f), e); continue
        for x in J:
            n = [re.sub(r"[^0-9A-Z]", "", str(t).upper().split("/")[-1]) for t in (x.get("numbers") or [])]; d[str(x.get("row")).lstrip("#")] = {"n": [t for t in n if t], "sure": bool(x.get("sure", True)), "note": str(x.get("note", "")), "form": str(x.get("form", ""))}
    return d
A, B = load("A"), load("B"); print("reader A panels", len(A), "| reader B panels", len(B))
def settle(rows, machine, layer_no):
    out = []
    for r in rows:
        a, b = A.get(r), B.get(r); m = [t for t in machine.get(r, "").split("/") if t]; L0 = layer_no.get(r, "")
        if a is None and b is None: out.append((r, L0, "not read yet", False)); continue
        if a is None or b is None:                                   # one reader only: it counts when it agrees with the second (machine) read
            o = a or b
            if len(o["n"]) == 1 and o["sure"] and o["n"][0] in m: out.append((r, o["n"][0], "one reader and the second read agree", True))
            else: out.append((r, L0, "one reader only: %s" % ("/".join(o["n"]) or "none"), False))
            continue
        na, nb = a["n"], b["n"]
        if na == nb and len(na) == 1: out.append((r, na[0], "two readers agree" + (" with the second read" if na[0] in m else ""), True))
        elif na == nb and len(na) == 0: out.append((r, "", "two readers: no number printed in the plot", True))
        elif len(na) == 1 and na[0] in m and (len(nb) != 1 or not b["sure"]): out.append((r, na[0], "one reader and the second read agree", True))
        elif len(nb) == 1 and nb[0] in m and (len(na) != 1 or not a["sure"]): out.append((r, nb[0], "one reader and the second read agree", True))
        elif len(na) == 1 and len(nb) == 1 and na[0] != nb[0] and (na[0] in m) != (nb[0] in m): out.append((r, na[0] if na[0] in m else nb[0], "one reader and the second read agree (the other reader saw %s)" % (nb[0] if na[0] in m else na[0]), True))
        else: out.append((r, L0, "readers differ: %s | %s" % ("/".join(na) or "none", "/".join(nb) or "none"), False))
    return pd.DataFrame(out, columns=["row", "final", "how", "settled"])
V = pd.read_csv(os.path.join(HERE, "out", "visual.csv"), dtype=str).fillna("")
X = settle(list(V.row), dict(zip(V.row, V.second_read + "/" + V.beside + "/" + V.first_read)), dict(zip(V.row, V.subdiv_no)))
print("rows looked at", len(X), "| settled", int(X.settled.sum()), "|", X.how.str.replace(r" \(.*| \|.*|: .*", "", regex=True).value_counts().to_dict())
X[X.settled].to_csv(os.path.join(HERE, "out", "visual_read.csv"), index=False); X[~X.settled].to_csv(os.path.join(HERE, "out", "visual_open.csv"), index=False)
if os.path.exists(os.path.join(HERE, "out", "visual2.csv")):
    V2 = pd.read_csv(os.path.join(HERE, "out", "visual2.csv"), dtype=str).fillna("")
    X2 = settle(list(V2.row), dict(zip(V2.row, V2.second_read)), {r: "" for r in V2.row})
    print("pieces looked at", len(X2), "| settled", int(X2.settled.sum()), "|", X2.how.str.replace(r" \(.*| \|.*|: .*", "", regex=True).value_counts().to_dict())
    X2.to_csv(os.path.join(HERE, "out", "visual_read2.csv"), index=False)
