import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd, math
from shapely.ops import unary_union
B = gpd.read_file(sys.argv[1], layer="plan").to_crs(32644)
for s in sys.argv[2:]:
    d = unary_union(list(B.geometry[B.survey == s])).buffer(0.01).buffer(-0.01).simplify(0.05)
    c = np.array(d.exterior.coords)[:-1]
    out = []
    for i in range(len(c)):
        a, b, e = c[i-1], c[i], c[(i+1) % len(c)]
        ang = math.degrees(math.acos(np.clip(np.dot(a-b, e-b) / np.linalg.norm(a-b) / np.linalg.norm(e-b), -1, 1)))
        if ang < 172: out.append("(%.1f, %.1f) %.0fdeg next %.1fm" % (b[0], b[1], ang, np.linalg.norm(e-b)))
    print(s, "area %.0f" % d.area); [print("   ", o) for o in out]
