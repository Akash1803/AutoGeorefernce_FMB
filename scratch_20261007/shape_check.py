"""Read-only proper check: every placed Paranur parcel against its own PDF drawing.
The PDF outline is fitted rigidly onto the parcel (vertex-to-boundary ICP from 36 start headings); every stretch of the
parcel's boundary that leaves the PDF line by more than TOL_M is reported with its length, largest offset and place.
Usage: python shape_check.py [disk|<plan.gpkg>]"""
import sys, os, math, json
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, shapely
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import unary_union
from autogeoref import ring_io, ring, fit, printed_numbers as pn, paths
TOL_M = 0.3
V = "35_04_221"
src = sys.argv[1] if len(sys.argv) > 1 and __name__ == "__main__" else "disk"

def parcels():
    out = {}
    if src == "disk":
        for L in ring_io.village_layers("Paranur"):
            n = L["name"]
            if "_parcels" not in n:
                continue
            s = n.split("_parcels")[0]
            if s in ("70A", "70B") and not n.endswith("_auto"):
                continue
            out[s] = (n, ring_io._read_source(L["source"]).to_crs(32644))
    else:
        P = gpd.read_file(src, layer="plan").to_crs(32644)
        for s, g in P.groupby("survey"):
            out[s] = (g.layer.iat[0], g)
    return out

def fit_pdf(sheet_outline, target):
    """Best rigid pose of the PDF outline onto the target polygon's boundary."""
    S = np.asarray(sheet_outline); T = target.exterior
    cs, ct = S.mean(0), np.asarray(target.centroid.coords[0])
    best = None
    for th0 in range(0, 360, 10):
        th, t = float(th0), ct - fit.transform_points([cs], th0, (0, 0))[0]
        for _ in range(25):
            P = fit.transform_points(S, th, t)
            Q = np.array([shapely.ops.nearest_points(T, Point(p))[0].coords[0] for p in P])
            th, t, rms, mx = fit.rigid_fit(S, Q)
        G = Polygon(fit.transform_points(S, th, t))
        iou = G.intersection(target).area / G.union(target).area if G.is_valid else 0
        if best is None or iou > best[0]:
            best = (iou, th, t)
    return best

def main():
    table = pn.village_table(V, str(paths.PROJECT / "Auto_Georef" / V / "printed_numbers.json"))
    rows, devs = [], []
    for s, (layer, g) in sorted(parcels().items(), key=lambda kv: paths.survey_sort_key(kv[0])):
        sh = ring.load_sheet(V, s, table)
        target = ring.largest(unary_union(g.geometry.values).buffer(0))
        iou, th, t = fit_pdf(sh.outline, target)
        pdf = Polygon(fit.transform_points(sh.outline, th, t))
        # boundary of the parcel sampled every 0.25 m: distance to the PDF line
        pts, _b = ring.sheets.sample_outline(list(target.exterior.coords)[:-1], 0.25)
        d = shapely.distance(shapely.points(pts), pdf.exterior)
        off = d > TOL_M
        # group consecutive off-line samples into stretches
        stretches, cur = [], []
        for k in range(len(pts)):
            if off[k]:
                cur.append(k)
            elif cur:
                stretches.append(cur); cur = []
        if cur:
            if stretches and stretches[0][0] == 0:
                stretches[0] = cur + stretches[0]
            else:
                stretches.append(cur)
        # PDF corners missing from the parcel
        pdf_pts = np.array(pdf.exterior.coords)[:-1]
        miss = shapely.distance(shapely.points(pdf_pts), target.exterior)
        for st in stretches:
            k = st[int(np.argmax(d[st]))]
            ll = gpd.GeoSeries([Point(pts[k])], crs=32644).to_crs(4326).iat[0]
            devs.append({"survey": s, "layer": layer, "length_m": round(0.25 * len(st), 1), "max_offset_m": round(float(d[st].max()), 2),
                         "lon": round(ll.x, 7), "lat": round(ll.y, 7)})
        rows.append({"survey": s, "layer": layer, "pdf_iou": round(iou, 4), "area_vs_pdf_%": round(100 * (target.area / pdf.area - 1), 2),
                     "boundary_off_pdf_m": round(0.25 * off.sum(), 1), "stretches_off_pdf": len(stretches),
                     "max_offset_m": round(float(d.max()), 2), "pdf_corners_missing": int((miss > TOL_M).sum())})
    R = pd.DataFrame(rows); D = pd.DataFrame(devs)
    pd.set_option("display.width", 220)
    print(R.to_string(index=False))
    print("\nlargest departures from the PDF line:")
    print(D.sort_values("max_offset_m", ascending=False).head(25).to_string(index=False))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "shape_check_%s.csv" % ("disk" if src == "disk" else "plan"))
    D.to_csv(out, index=False); R.to_csv(out.replace(".csv", "_summary.csv"), index=False)


if __name__ == "__main__":
    main()
