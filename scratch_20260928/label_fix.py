"""Plot-number corrections from the PDF cross-check of 2026-09-28 -> out/label_fix.json {file|layer: {fid: {...}}}
Every change was read on the sheet (character outlines) and looked at on a crop of the PDF."""
import warnings, sys, os, json; warnings.filterwarnings("ignore")
sys.path.insert(0, r"D:\code\FMB_to_GeoJSON")
import pandas as pd, geopandas as gpd
from autogeoref import visible
M = pd.read_csv("out/decide4.csv", dtype=str).fillna("").set_index("rid")
SET = {  # rid: new plot number
    "50": "4", "81": "2G2", "82": "2G1", "159": "1C1A2", "504": "2B2A",
    "255": "6A", "257": "1H1A", "258": "1E", "261": "73", "262": "1D2A1A", "272": "1H1A1", "278": "1J1A1A1A1", "279": "1A1", "289": "1I1A1A", "321": "63", "326": "101",
    "343": "124", "344": "1J1A2", "346": "107", "348": "1J1B", "361": "1H2", "365": "66", "369": "52", "373": "1D2B", "375": "1F2", "377": "53B", "384": "70A", "387": "46A",
    "388": "108", "391": "110", "392": "114", "399": "99", "400": "106", "401": "1I2", "403": "54A1", "404": "55A1", "405": "54B", "406": "109", "407": "55B", "408": "55A2",
    "409": "54A2", "411": "1G1A", "414": "113", "416": "131", "417": "1G1B",
    "548": "2B1", "554": "2C3D", "563": "2C2D", "565": "2C2E", "566": "5", "567": "2C2H", "570": "2C2F", "572": "2C2G", "606": "2A1D",
    "610": "1A1E", "622": "1C1D", "623": "1C1E", "631": "1B1D", "635": "1C1F", "638": "1A1D"}
CLEAR = ["259", "270", "271", "285", "286", "294", "332", "364"]           # 526A: the '6' was a digit of a dimension, no number is printed
WHY = {"81": "printed beside the plot among overlapping text; 2G2 is the label nearest this plot", "82": "printed beside the plot among overlapping text; 2G1 is the label nearest this plot",
       "391": "printed as 526/110", "566": "printed 5, overlapped by the dimension 10.4; it lies between plots 4 and 6", "50": "printed 4, overlapped by the dimension (92.2)",
       "289": "printed 1I1A1A (the letter I was not read before)", "401": "printed 1I2 (the letter I was not read before)", "159": "printed 1C1A2 (the last 2 is overlapped by a dimension)",
       "278": "printed 1J1A1A1A1 across this strip", "392": "printed 114, overlapped by the dimension 14.0"}
for k in SET: assert k in M.index, k
# no number may repeat inside a sheet afterwards
new = M.layer_plot.copy()
for k, v in SET.items(): new[k] = v
for k in CLEAR: assert M.layer_plot[k] == "6" and M.survey[k] == "526A", k; new[k] = ""
chk = pd.DataFrame({"v": M.village_code, "s": M.survey, "n": new, "single": M.single}); chk = chk[(chk.n != "") & (chk.single == "False")]
dup = chk[chk.duplicated(["v", "s", "n"], keep=False)]; assert dup.empty, dup
plan = {}
for V in ["35_04_074", "35_04_077"]:
    sh = visible.shown(V)
    for s, (f, l) in sh.items():
        rows = M[(M.village_code == V) & (M.survey == s)]
        if not any(r in SET or r in CLEAR for r in rows.index): continue
        L = gpd.read_file(f, layer=l, fid_as_index=True); pos = gpd.read_file(f, layer=l)            # same row order: real FIDs and the row numbers the ledger used
        for rid, r in rows.iterrows():
            if rid not in SET and rid not in CLEAR: continue
            fid = int(L.index[int(r.fid)]); old = pos.plot_no.iloc[int(r.fid)]; old = "" if old is None or str(old) in ("nan", "None") else str(old)
            assert old == r.layer_plot, (s, rid, old, r.layer_plot); assert abs(L.geometry.iloc[int(r.fid)].area - float(r.area)) < 1.0 or L.crs.to_epsg() != 32644
            if rid in SET:
                src = "checked on sheet 2026-09-28: " + WHY.get(rid, "printed %s" % SET[rid]) + (" (was %s)" % old if old else "")
                plan.setdefault(f"{f}|{l}", {})[str(fid)] = {"plot_no": SET[rid], "numbered": True, "review_needed": rid in ("81", "82"), "label_source": src, "old": old}
            else:
                plan.setdefault(f"{f}|{l}", {})[str(fid)] = {"plot_no": None, "numbered": False, "review_needed": True, "label_source": "checked on sheet 2026-09-28: no number printed (the 6 read earlier is a digit of a dimension)", "old": old}
json.dump(plan, open("out/label_fix.json", "w"), indent=1)
print({os.path.basename(k.split("|")[0]): len(v) for k, v in plan.items()}, "| total", sum(len(v) for v in plan.values()))
