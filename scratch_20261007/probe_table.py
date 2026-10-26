"""Probe: printed neighbour numbers on every sheet of a village, and two-way agreement."""
import sys, os, re, json, glob, math
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import pdfplumber
import fmb_to_geojson as C
from autogeoref import neighbours
from probe_black import black_glyphs, templates
from probe_runs import glyphs_of, runs_of, read
V = sys.argv[1]; PR = r"D:\Projects\Tambaram_Chengalpattu"
T = templates()
def scale_of(s):
    f = os.path.join(PR, "FMB_Vector", V, s + "_parcels.geojson")
    sc = json.load(open(f, encoding="utf-8"))["features"][0]["properties"].get("sheet_scale", "1:500")
    return int(str(sc).split(":")[-1])
table = {}
for pdf in sorted(glob.glob(os.path.join(PR, "FMB_Sketches", V, "*.pdf"))):
    s = os.path.splitext(os.path.basename(pdf))[0]
    if not os.path.exists(os.path.join(PR, "FMB_Vector", V, s + "_parcels.geojson")): continue
    page = pdfplumber.open(pdf).pages[0]; H = float(page.height); W = float(page.width)
    M = C.PT_MM * scale_of(s) / 1000.0
    gl = glyphs_of(black_glyphs(page)); out = []
    for r in runs_of(gl):
        x0 = min(gl[k]["x0"] for k in r); x1 = max(gl[k]["x1"] for k in r); y0 = min(gl[k]["y0"] for k in r); y1 = max(gl[k]["y1"] for k in r)
        cy = (y0 + y1) / 2
        if cy < 0.04 * H or cy > 0.972 * H: continue
        txt, worst = read(gl, r, T)
        if worst > 40 or not re.fullmatch(r"\d+[A-Z]?(/[A-Z0-9]+)*", txt): out.append(("?" + txt, worst)); continue
        out.append((txt, worst))
    table[s] = out
nb = {s: {neighbours.normalise(t) for t, w in v if not t.startswith("?")} - {None} for s, v in table.items()}
for s, v in table.items():
    print("%-5s prints %s" % (s, ", ".join("%s(%d)" % (t, w) for t, w in v)))
print("\n-- two-way agreement among sheets we have --")
for s in table:
    both = sorted(n for n in nb[s] if n in nb and s.upper() in nb[n])
    one = sorted(n for n in nb[s] if n in nb and s.upper() not in nb[n])
    print("%-5s agreed: %-40s one-sided: %s" % (s, ",".join(both), ",".join(one)))
