import numpy as np, rasterio, cv2, geopandas as gpd, glob, os, math, warnings
warnings.filterwarnings("ignore")
from shapely.geometry import LineString, Point
from shapely.strtree import STRtree
SP=r"C:\Users\FAI-Akash\AppData\Local\Temp\claude\d--Akash-Python\be3d1523-acf3-4fef-9ca8-783aeaecb2d3\scratchpad"
tif=os.path.join(SP,"spike_gcp","kizhi_sat_z20.tif")
with rasterio.open(tif) as ds:
    img=ds.read(); T=ds.transform
gray=cv2.cvtColor(np.moveaxis(img,0,-1),cv2.COLOR_RGB2GRAY)
lsd=cv2.createLineSegmentDetector(cv2.LSD_REFINE_STD)
lines=lsd.detect(gray)[0].reshape(-1,4)
def px2map(x,y): return (T.c+x*T.a, T.f+y*T.e)
segs=[]; segarr=[]
for x1,y1,x2,y2 in lines:
    p1=px2map(x1,y1); p2=px2map(x2,y2)
    L=math.hypot(p2[0]-p1[0],p2[1]-p1[1])
    if L<3: continue
    segs.append(LineString([p1,p2])); segarr.append((np.array(p1),np.array(p2)))
tree=STRtree(segs)
print("segments>=3m:",len(segs))

def outline(gdf):
    g=gdf.geometry.union_all().buffer(0)
    if g.geom_type=="MultiPolygon": g=max(g.geoms,key=lambda p:p.area)
    return g.exterior

def angdiff(a,b):
    d=abs((a-b)%math.pi); return min(d,math.pi-d)

def edge_matches(c, v, band=2.5, dth=math.radians(20)):
    """segments matched to edge c->v: direction within dth, perpendicular dist<=band over an overlap>=2 m."""
    c=np.array(c); v=np.array(v); e=v-c; L=np.linalg.norm(e); u=e/L; n=np.array([-u[1],u[0]])
    th=math.atan2(u[1],u[0])%math.pi
    edge=LineString([tuple(c),tuple(v)])
    out=[]
    for j in tree.query(edge.buffer(band)):
        p1,p2=segarr[j]
        sth=math.atan2(*(p2-p1)[::-1])%math.pi
        if angdiff(sth,th)>dth: continue
        a1=(p1-c)@u; a2=(p2-c)@u; lo=max(0,min(a1,a2)); hi=min(L,max(a1,a2)); ov=hi-lo
        if ov<1.0: continue
        pa=p1+(p2-p1)*((lo-a1)/(a2-a1)) if a2!=a1 else p1; pb=p1+(p2-p1)*((hi-a1)/(a2-a1)) if a2!=a1 else p2; d1=abs((pa-c)@n); d2=abs((pb-c)@n)
        if max(d1,d2)>band: continue
        out.append((j,ov,lo,(d1+d2)/2))
    return out

def line_intersect(s1,s2):
    p,r=s1[0],s1[1]-s1[0]; q,s=s2[0],s2[1]-s2[0]
    den=r[0]*s[1]-r[1]*s[0]
    if abs(den)<1e-9: return None
    t=((q-p)[0]*s[1]-(q-p)[1]*s[0])/den
    return p+t*r

def coverage(c,v,band=1.5,dth=math.radians(12)):
    """share of edge length within band of a same-direction segment (the spike metric), per edge."""
    c=np.array(c); v=np.array(v); L=np.linalg.norm(v-c); n=max(int(L),2)
    th=math.atan2(*(v-c)[::-1])%math.pi; hit=0
    for i in range(n):
        p=Point(c+(v-c)*(i+0.5)/n)
        for j in tree.query(p.buffer(band)):
            p1,p2=segarr[j]; sth=math.atan2(*(p2-p1)[::-1])%math.pi
            if angdiff(sth,th)<=dth and segs[j].distance(p)<=band: hit+=1; break
    return hit/n

bk=r"D:\Projects\Tambaram_Chengalpattu\_logs\backup_georef_20260917"
for f in sorted(glob.glob(os.path.join(bk,"*.gpkg"))):
    sv=os.path.basename(f).split("_")[0]
    gdf=gpd.read_file(f)
    if gdf.crs is None or gdf.crs.to_epsg()!=32644: gdf=gdf.to_crs(32644)
    ring=outline(gdf); ring=ring.simplify(0.05)
    allpts=[]
    for gm in gdf.geometry:
        for pg in (gm.geoms if gm.geom_type=="MultiPolygon" else [gm]): allpts.append(list(pg.exterior.simplify(0.05).coords)[:-1])
    counts={1.0:set(),2.0:set(),3.0:set()}
    corners=0; both=0; edge_cov=[]
    rings=[list(ring.coords)[:-1]]+allpts
    for ri,pts in enumerate(rings):
     N=len(pts)
     for i in range(N):
         pv=np.array(pts[i-1]); c=np.array(pts[i]); nx=np.array(pts[(i+1)%N])
         a=math.atan2(*(c-pv)[::-1]); b=math.atan2(*(nx-c)[::-1]); turn=abs((b-a+math.pi)%(2*math.pi)-math.pi)
         if math.degrees(turn)<20: continue
         corners+=1
         m1=edge_matches(c,pv); m2=edge_matches(c,nx)
         if not m1 or not m2: continue
         both+=1
         # try every pair among the top candidates (generous): best overlap, nearest to corner
         cand1=sorted(m1,key=lambda x:-x[1])[:3]+sorted(m1,key=lambda x:x[2])[:3]
         cand2=sorted(m2,key=lambda x:-x[1])[:3]+sorted(m2,key=lambda x:x[2])[:3]
         bestd=1e9
         for s1 in cand1:
             for s2 in cand2:
                 X=line_intersect(segarr[s1[0]],segarr[s2[0]])
                 if X is None: continue
                 bestd=min(bestd,np.linalg.norm(X-c))
         for t in counts:
             if bestd<=t: counts[t].add((ri,i))
    # per-edge coverage for diagnostics
    for i in range(N):
        c=np.array(pts[i]); nx=np.array(pts[(i+1)%N]); L=np.linalg.norm(nx-c)
        if L>=5: edge_cov.append((L,coverage(c,nx)))
    pts=list(ring.coords)[:-1]
    def span(idx):
        if len(idx)<2: return 0.0
        P=np.array([rings[r][i] for r,i in idx]); return max(np.linalg.norm(P[a]-P[b]) for a in range(len(P)) for b in range(a+1,len(P)))
    uniq={t:len({(round(rings[r][i][0],1),round(rings[r][i][1],1)) for r,i in counts[t]}) for t in counts}
    print("   unique GCP corners @1m/@2m/@3m:",uniq[1.0],uniq[2.0],uniq[3.0])
    ec=" ".join(f"{L:.0f}m:{cv*100:.0f}%" for L,cv in edge_cov)
    print(f"{sv:4s} corners(turn>=20)={corners:2d} both-edges-matched={both:2d} | GCPs @1m={len(counts[1.0])} (span {span(counts[1.0]):.0f} m) @2m={len(counts[2.0])} @3m={len(counts[3.0])} | edges>=5m coverage: {ec}",flush=True)
