"""Dry run: place 70A again, (a) with the two-way rule, (b) also trusting the neighbours 70A's own sheet prints
(72 prints nothing back). Reports overlaps / gaps of each against everything else. Nothing is written."""
import sys
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, pandas as pd, shapely
from shapely.ops import unary_union
from shapely.geometry import Polygon
from autogeoref import ring, ring_io, ring_topo, printed_numbers as pn, paths
V = "35_04_221"
seeds, sources, plots = ring_io.load_seeds("Paranur", include_red=True)
for s in ("70A",):
    seeds.pop(s, None)
plots = plots[plots.survey_no != "70A"]
table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
print("70A prints:", sorted(pn.printed_sets(table)["70A"]), "| two-way:", sorted(pn.agreed_sets(table)["70A"]))
orig_named = ring.Context.named
for mode in ("two-way", "own sheet"):
    if mode == "own sheet":
        def named(self, a, b):
            r = orig_named(self, a, b)
            if r == "one-sided" and str(b).upper() in self.printed.get(str(a).upper(), set()):
                return "agreed"
            return r
        ring.Context.named = named
    dec, sheets_, ctx = ring.place_ring(V, seeds, table, ring=["70A"], in_order=True, log=lambda *a: None)
    d = dec[0]
    if d.status != "placed":
        print(mode, "-> not placed:", d.reason); continue
    G = ring.ground_polygon(sheets_["70A"], d.theta, d.t)
    U = {s: g for s, g in seeds.items()}
    ov = {s: round(G.intersection(g).area, 1) for s, g in U.items() if G.intersection(g).area > 0.5}
    allu = unary_union(list(U.values()) + [G])
    holes = [round(Polygon(r).difference(allu).area, 1) for p in getattr(allu, "geoms", [allu]) for r in p.interiors]
    print("%-9s support %s | overlaps %s | gaps %s | 70A to 72 %.2f m" % (mode, d.metrics.get("support_by"), ov, [h for h in holes if h > 1], G.distance(U["72"])))
