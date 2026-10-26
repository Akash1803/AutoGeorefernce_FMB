"""Decision table from the position-matched ledger. -> out/decide4.csv"""
import pandas as pd
M = pd.read_csv("out/ledger4.csv", dtype=str).fillna(""); M["area"] = M.area.astype(float); M["w_iou"] = M.w_iou.astype(float)
claim = M[M.single == "False"].groupby(["village_code", "survey", "sheet_pid"]).w_iou.transform("max"); M["loser"] = False
M.loc[claim.index, "loser"] = (M.loc[claim.index, "w_iou"] < claim) & (M.loc[claim.index, "w_iou"] < 0.3)
used = M[M.layer_plot != ""].groupby(["village_code", "survey"]).layer_plot.apply(list).to_dict()
st, prop = [], []
for r in M.itertuples():
    L = r.layer_plot; c = [x for x in r.labels.split("/") if x]; noise6 = r.survey == "526A" and L == "6"
    strong = [x for x in c if len(x) >= 2 and x != "00"]
    if r.single == "True": st.append("whole survey"); prop.append(""); continue
    if r.loser: st.append("sliver"); prop.append(""); continue
    if L and not noise6 and L in c: st.append("confirmed"); prop.append(L)
    elif L and not noise6 and c: st.append("conflict"); prop.append("/".join(c))
    elif L and not noise6: st.append("unverified"); prop.append(L)
    elif len(strong) == 1: st.append("proposed"); prop.append(strong[0])
    elif len(strong) > 1: st.append("several"); prop.append("/".join(strong))
    else: st.append("none"); prop.append("/".join(c))
M["status"], M["proposed"] = st, prop; M.to_csv("out/decide4.csv", index=False)
print(M.status.value_counts().to_dict()); print("in buffer:", M[M.in_buffer == "True"].status.value_counts().to_dict())
for s_ in ["conflict", "unverified", "proposed", "several", "sliver"]:
    x = M[M.status == s_]; print("\n==", s_, len(x))
    for (v, s), g in x.groupby(["village_code", "survey"], sort=False):
        print(v[-3:], s, "; ".join(f"#{r.rid} {r.layer_plot or '?'}->{r.proposed}({r.area:.0f}m2,{r.w_iou:.2f}){' B' if r.in_buffer == 'True' else ''}" for r in g.itertuples()))
x = M[M.status == "none"]; print("\n== none", len(x))
for (v, s), g in x.groupby(["village_code", "survey"], sort=False):
    print(v[-3:], s, len(g), "; ".join(f"#{r.rid}({r.area:.0f}m2{',' + r.proposed if r.proposed else ''}){' B' if r.in_buffer == 'True' else ''}" for r in g.itertuples()))
