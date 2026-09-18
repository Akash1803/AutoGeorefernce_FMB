import numpy as np, rasterio, cv2, geopandas as gpd, glob, os, math, warnings, time
warnings.filterwarnings("ignore")
SP=r"C:\Users\FAI-Akash\AppData\Local\Temp\claude\d--Akash-Python\be3d1523-acf3-4fef-9ca8-783aeaecb2d3\scratchpad"
tif=os.path.join(SP,"spike_gcp","kizhi_sat_z20.tif")
with rasterio.open(tif) as ds:
    img=ds.read(); T=ds.transform; b=ds.bounds
gray=cv2.cvtColor(np.moveaxis(img,0,-1),cv2.COLOR_RGB2GRAY)
lsd=cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD)
lines=lsd.detect(gray)[0].reshape(-1,4)
def px2map(x,y): return (T.c+x*T.a, T.f+y*T.e)
segs=[]
for x1,y1,x2,y2 in lines:
    p1=px2map(x1,y1); p2=px2map(x2,y2)
    L=math.hypot(p2[0]-p1[0],p2[1]-p1[1])
    if L<3: continue
    segs.append((p1,p2,math.atan2(p2[1]-p1[1],p2[0]-p1[0])%math.pi))
print("segments>=3m",len(segs),"window m2 %.0f"%((b.right-b.left)*(b.top-b.bottom)),flush=True)
CELL=0.5; NB=12
W=int(math.ceil((b.right-b.left)/CELL)); H=int(math.ceil((b.top-b.bottom)/CELL))
dts=np.zeros((NB,H,W),np.float32)
for k in range(NB):
    mask=np.full((H,W),255,np.uint8)
    for p1,p2,th in segs:
        if int(th/(math.pi/NB))%NB!=k: continue
        c1=((p1[0]-b.left)/CELL,(b.top-p1[1])/CELL); c2=((p2[0]-b.left)/CELL,(b.top-p2[1])/CELL)
        cv2.line(mask,(int(round(c1[0])),int(round(c1[1]))),(int(round(c2[0])),int(round(c2[1]))),0,1)
    dts[k]=cv2.distanceTransform(mask,cv2.DIST_L2,5)*CELL
binc=(np.arange(NB)+0.5)*(math.pi/NB)
SIG=0.75; BAND=2.5
def score(xy,th):
    cx=((xy[:,0]-b.left)/CELL).astype(int); cy=((b.top-xy[:,1])/CELL).astype(int)
    ok=(cx>=0)&(cx<W)&(cy>=0)&(cy<H); cx=np.clip(cx,0,W-1); cy=np.clip(cy,0,H-1)
    d=dts[:,cy,cx]; dth=np.abs(binc[:,None]-th[None,:]); dth=np.minimum(dth,math.pi-dth)
    s=np.exp(-d*d/(2*SIG*SIG))*np.cos(dth)**2; s[d>BAND]=0
    best=s.max(axis=0); best[~ok]=0
    return float(best.mean())
def outline(gdf):
    g=gdf.geometry.union_all().buffer(0)
    if g.geom_type=="MultiPolygon": g=max(g.geoms,key=lambda p:p.area)
    return g.exterior
def samples(ring):
    n=max(int(ring.length),4); L=ring.length
    P=np.array([ring.interpolate(i*L/n).coords[0] for i in range(n)])
    Q1=np.array([ring.interpolate(((i+0.5)*L/n)%L).coords[0] for i in range(n)])
    Q0=np.array([ring.interpolate(((i-0.5)*L/n)%L).coords[0] for i in range(n)])
    return P,(np.arctan2(Q1[:,1]-Q0[:,1],Q1[:,0]-Q0[:,0])%math.pi)
def posed(P,th,c,dth_deg,dx,dy):
    a=math.radians(dth_deg); R=np.array([[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]])
    return (P-c)@R.T+c+np.array([dx,dy]), (th+a)%math.pi
bk=r"D:\Projects\Tambaram_Chengalpattu\_logs\backup_georef_20260917"
ROT=range(-10,11,1); SH=range(-20,21,1)
print("svy   len_m  true  best_wrong (dx,dy,dth)  n_poses>=0.95true  ridge_extent_m(along/across)  rank_of_true",flush=True)
for f in sorted(glob.glob(os.path.join(bk,"*.gpkg"))):
    sv=os.path.basename(f).split("_")[0]
    gdf=gpd.read_file(f)
    if gdf.crs is None or gdf.crs.to_epsg()!=32644: gdf=gdf.to_crs(32644)
    ring=outline(gdf); c=np.array(ring.centroid.coords[0]); P,th=samples(ring)
    s_true=score(P,th)
    res=[]
    for r in ROT:
        for dx in SH:
            for dy in SH:
                s=score(*posed(P,th,c,r,dx,dy)); res.append((s,dx,dy,r))
    res=np.array(res)
    wrong=res[(np.hypot(res[:,1],res[:,2])>=3)|(np.abs(res[:,3])>=2)]
    bw=wrong[wrong[:,0].argmax()]
    near=res[res[:,0]>=0.95*s_true]
    # ridge extent: spread of near-max poses along principal axis
    if len(near)>1:
        xy=near[:,1:3]; xy=xy-xy.mean(0); u,sv_,vt=np.linalg.svd(xy,full_matrices=False)
        proj=xy@vt[0]; proj2=xy@vt[1]; ext=(proj.max()-proj.min(), proj2.max()-proj2.min()); bearing=math.degrees(math.atan2(vt[0][0],vt[0][1]))%180
    else: ext=(0,0); bearing=float('nan')
    rank=int((res[:,0]>s_true).sum())
    print(f"{sv:5s} {ring.length:6.0f} {s_true:5.2f}  {bw[0]:5.2f} ({bw[1]:+.0f},{bw[2]:+.0f},{bw[3]:+.0f})   {len(near):5d}   {ext[0]:5.1f}/{ext[1]:5.1f} bearing {bearing:5.1f}   {rank}",flush=True)
