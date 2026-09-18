import numpy as np, rasterio, cv2, geopandas as gpd, glob, os, math
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree
from shapely.affinity import rotate, translate
SP=r"C:\Users\FAI-Akash\AppData\Local\Temp\claude\d--Akash-Python\be3d1523-acf3-4fef-9ca8-783aeaecb2d3\scratchpad"
tif=os.path.join(SP,"spike_gcp","kizhi_sat_z20.tif")
with rasterio.open(tif) as ds:
    img=ds.read(); T=ds.transform
gray=cv2.cvtColor(np.moveaxis(img,0,-1),cv2.COLOR_RGB2GRAY)
lsd=cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD)
lines=lsd.detect(gray)[0].reshape(-1,4)
def px2map(x,y): return (T.c+x*T.a, T.f+y*T.e)
segs=[];dirs=[]
for x1,y1,x2,y2 in lines:
    p1=px2map(x1,y1); p2=px2map(x2,y2)
    L=math.hypot(p2[0]-p1[0],p2[1]-p1[1])
    if L<3: continue
    segs.append(LineString([p1,p2])); dirs.append(math.atan2(p2[1]-p1[1],p2[0]-p1[0])%math.pi)
dirs=np.array(dirs); tree=STRtree(segs)
print("segments >=3 m:",len(segs), "total length m: %.0f"%sum(s.length for s in segs))
def outline(gdf):
    g=gdf.geometry.union_all() if hasattr(gdf.geometry,"union_all") else gdf.geometry.unary_union
    g=g.buffer(0)
    if g.geom_type=="MultiPolygon": g=max(g.geoms,key=lambda p:p.area)
    return g.exterior
def score(ring, sigma=0.75, band=2.5, spike=False):
    n=max(int(ring.length),4); pts=[ring.interpolate(i*ring.length/n) for i in range(n)]
    # local direction of ring at sample
    sc=[]
    for i,p in enumerate(pts):
        q=ring.interpolate(((i+0.5)*ring.length/n)%ring.length); q0=ring.interpolate(((i-0.5)*ring.length/n)%ring.length)
        th=math.atan2(q.y-q0.y,q.x-q0.x)%math.pi
        idx=tree.query(p.buffer(band))
        best=0.0
        for j in idx:
            d=segs[j].distance(p)
            if d>band: continue
            dth=abs(dirs[j]-th); dth=min(dth,math.pi-dth)
            if spike:
                s=1.0 if (d<=1.5 and dth<=math.radians(12)) else 0.0
            else:
                s=math.exp(-d*d/(2*sigma*sigma))*math.cos(dth)**2
            best=max(best,s)
        sc.append(best)
    return float(np.mean(sc))
bk=os.path.join(r"D:\Projects\Tambaram_Chengalpattu\_logs\backup_georef_20260917")
rows=[]
for f in sorted(glob.glob(os.path.join(bk,"*.gpkg"))):
    sv=os.path.basename(f).split("_")[0]
    gdf=gpd.read_file(f)
    if gdf.crs is None or gdf.crs.to_epsg()!=32644: gdf=gdf.to_crs(32644)
    ring=outline(gdf); c=ring.centroid
    s_true=score(ring); s_spike=score(ring,spike=True)
    # wrong poses: shifts on a 2 m grid within +-20 m, rotations 0,+-3,+-6 deg; keep those >=3 m from truth
    best_wrong=(0,None)
    for dth in (0,-3,3,-6,6):
        r=rotate(ring,dth,origin=c) if dth else ring
        for dx in range(-20,21,2):
            for dy in range(-20,21,2):
                if math.hypot(dx,dy)<3 and dth==0: continue
                s=score(translate(r,dx,dy))
                if s>best_wrong[0]: best_wrong=(s,(dx,dy,dth))
    rows.append((sv,ring.length,s_true,s_spike,best_wrong))
    print(f"{sv:5s} len={ring.length:6.1f} m  spec-score@manual={s_true:.2f}  spike-share@manual={s_spike:.2f}  best WRONG pose score={best_wrong[0]:.2f} at (dx,dy,dth)={best_wrong[1]}", flush=True)
