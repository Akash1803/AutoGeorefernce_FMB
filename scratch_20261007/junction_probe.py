"""Read-only: re-seat 72 and 70A on the 69/70A/72/73/74 junction the four PDFs agree on, and measure what is left."""
import sys; sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, shapely
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import unary_union
from autogeoref import ring_io, ring, fit, printed_numbers as pn, paths
from shape_check import fit_pdf
V = "35_04_221"
table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
L = {l["name"]: l["source"] for l in ring_io.village_layers("Paranur")}
def get(n): return ring.largest(unary_union(list(ring_io._read_source(L[n]).to_crs(32644).geometry)).buffer(0))
g69, g267, g68B, g73, g72, g70A = (get(n) for n in ["69_parcels_modified", "267_parcels_modified", "68B_parcels_modified", "73_parcels_auto", "72_parcels_auto", "70A_parcels_auto"])
J = np.array([389698.1, 1407893.3]); B69 = np.array([389740.0, 1408008.2]); C73 = np.array([389604.7, 1407585.4])
def snap(p, g):
    c = np.array(g.exterior.coords)[:-1]; return c[np.argmin(np.linalg.norm(c - p, axis=1))]
J, B69, C73 = snap(J, g69), snap(B69, g69), snap(C73, g73)
print("J", J, "69.B", B69, "73.C", C73)
def seat(s, cur, pairs):
    sh = ring.load_sheet(V, s, table); iou, th, t = fit_pdf(sh.outline, cur)
    P = fit.transform_points(sh.outline, th, t)
    idx = [int(np.argmin(np.linalg.norm(P - np.array(g), axis=1))) for g, _ in pairs]
    S = np.array(sh.outline)[idx]; Q = np.array([q for _, q in pairs])
    th2, t2, rms, mx = fit.rigid_fit(S, Q)
    G = Polygon(fit.transform_points(sh.outline, th2, t2))
    print("%s: corners %s seated, rms %.2f m, shift from now %.2f m, turn %.2f deg" % (s, idx, rms, np.linalg.norm(np.array(G.centroid.coords[0]) - np.array(Polygon(P).centroid.coords[0])), th2 - th))
    return G
# 72: A (now at 389700.5,1407900.0) -> J ; G (now 389606.9,1407592.8) -> 73.C
G72 = seat("72", g72, [((389700.5, 1407900.0), J), ((389606.9, 1407592.8), C73)])
# 70A rigid: C (SW) -> J ; D (NW, along 69) -> 69.B
sh = ring.load_sheet(V, "70A", table); c70 = np.array(g70A.exterior.coords)[:-1]
cC = c70[np.argmin(np.linalg.norm(c70 - J, axis=1))]; cD = c70[np.argmin(np.linalg.norm(c70 - B69, axis=1) + 0 * c70[:, 0])]
# D of the rigid 70A is its western top corner: the corner farthest north on the 69 side
cD = max(c70, key=lambda p: p[1] - 3 * abs(LineString([J, B69]).distance(Point(p))))
G70 = seat("70A", g70A, [(tuple(cC), J), (tuple(cD), B69)])
line267 = g267.exterior
for nm, G in (("72", G72), ("70A", G70)):
    c = np.array(G.exterior.coords)[:-1]
    print(nm, "corners -> distance to 267 line:", [(tuple(np.round(p, 1)), round(line267.distance(Point(p)), 2)) for p in c])
    print(nm, "overlaps: 69 %.1f 267 %.1f 68B %.1f 73 %.1f | area %.0f" % (G.intersection(g69).area, G.intersection(g267).area, G.intersection(g68B).area, G.intersection(g73).area, G.area))
print("70A x 72 overlap %.1f" % G70.intersection(G72).area)
# the space left for 70A between 69, 267, 68B and the re-seated 72
U = unary_union([g69, g267, g68B, G72, g73])
box = Polygon([J + [-5, -5], B69 + [-5, 15], B69 + [40, 15], J + [60, -5]])
hole = box.difference(U)
pieces = sorted(getattr(hole, "geoms", [hole]), key=lambda q: -q.area)
sp = [q for q in pieces if q.distance(Point(J)) < 1 and q.distance(Point(B69)) < 1]
print("space between 69 / 267 / 68B / re-seated 72: %.0f m2 (70A sheet %.0f m2, register 2550)" % (sp[0].area if sp else -1, Polygon(sh.outline).area))
gpd.GeoDataFrame({"name": ["72_seated", "70A_seated"]}, geometry=[G72, G70], crs=32644).to_file("out/junction_probe.gpkg", driver="GPKG")
