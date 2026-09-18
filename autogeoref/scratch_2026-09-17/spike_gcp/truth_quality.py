import json, glob, os, math, warnings
warnings.filterwarnings("ignore")
import numpy as np, geopandas as gpd
from shapely.geometry import shape
from shapely.ops import unary_union
V = "D:/Projects/Tambaram_Chengalpattu/FMB_Vector/35_04_077/"
B = "D:/Projects/Tambaram_Chengalpattu/_logs/backup_georef_20260917/"
def fit(P,Q,scale):
    P=np.asarray(P,float);Q=np.asarray(Q,float);pc,qc=P.mean(0),Q.mean(0);A,Bm=P-pc,Q-qc
    H=A.T@Bm;U,S,Vt=np.linalg.svd(H);d=np.sign(np.linalg.det(Vt.T@U));D=np.diag([1,d]);R=Vt.T@D@U.T
    s=(S@np.diag([1,d])).sum()/(A**2).sum() if scale else 1.0
    t=qc-s*(R@pc);res=Q-(s*(R@P.T).T+t)
    return math.degrees(math.atan2(R[1,0],R[0,0])),s,t,math.sqrt((res**2).sum(1).mean()),np.sqrt((res**2).sum(1)).max(),res
for path in sorted(glob.glob(B+"*.gpkg")):
    s=os.path.basename(path).split("_")[0]
    g=gpd.read_file(path)
    if g.crs and g.crs.to_epsg()!=32644: g=g.to_crs(32644)
    d=json.load(open(V+f"{s}_parcels.geojson"))
    sheet={f['properties']['poly_id']:shape(f['geometry']) for f in d['features']}
    key='poly_id' if 'poly_id' in g.columns else None
    pairs=[]
    for pid,geom in sheet.items():
        row=g[g[key]==pid] if key else g.iloc[[pid-1]]
        if len(row)!=1: continue
        pg=row.geometry.iloc[0]
        if pg.geom_type=='MultiPolygon': pg=max(pg.geoms,key=lambda x:x.area)
        a=list(geom.exterior.coords)[:-1]; b=list(pg.exterior.coords)[:-1]
        if len(a)!=len(b): continue
        pairs+= [(pid,x,y) for x,y in zip(a,b)]
    if len(pairs)<3: print(s,"no pairs"); continue
    P=[p[1] for p in pairs];Q=[p[2] for p in pairs]
    th_r,_,t_r,rms_r,mx_r,res=fit(P,Q,False); th_s,sc,_,rms_s,mx_s,_=fit(P,Q,True)
    ext=unary_union(list(sheet.values())); minx,miny,maxx,maxy=ext.bounds; L=math.hypot(maxx-minx,maxy-miny)
    th_jv=[fit([p[1] for j,p in enumerate(pairs) if j!=i],[p[2] for j,p in enumerate(pairs) if j!=i],False)[0] for i in range(len(pairs))]
    polys=sorted({p[0] for p in pairs})
    th_jp=[fit([p[1] for p in pairs if p[0]!=k],[p[2] for p in pairs if p[0]!=k],False)[0] for k in polys] if len(polys)>1 else [th_r]
    man_c=unary_union(list(g.geometry)).centroid
    R=np.array([[math.cos(math.radians(th_r)),-math.sin(math.radians(th_r))],[math.sin(math.radians(th_r)),math.cos(math.radians(th_r))]])
    sc_c=np.array(ext.centroid.coords[0]); placed_c=R@sc_c+t_r
    dc=math.hypot(placed_c[0]-man_c.x,placed_c[1]-man_c.y)
    print(f"{s:4s} polys={len(sheet):2d} pairs={len(pairs):3d} extent={L:6.1f} m | rigid: heading {th_r:7.2f} rms {rms_r:.2f} max {mx_r:.2f} | similarity: scale {sc:.4f} heading {th_s:7.2f} rms {rms_s:.2f} | "
          f"jackknife heading spread: vertex {max(th_jv)-min(th_jv):.2f} deg, polygon {max(th_jp)-min(th_jp):.2f} deg | manual-centroid vs rigid-centroid {dc:.2f} m")
