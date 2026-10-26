import shapely
from shapely import make_valid
from shapely.geometry import Polygon
def valid4326(g):
    """a valid single Polygon for EPSG:4326 storage (micro self-intersections from reprojection/rounding removed)"""
    if g.is_valid and g.geom_type == "Polygon": return g
    for f in (lambda x: shapely.set_precision(x, 1e-9), lambda x: x.buffer(0), lambda x: make_valid(x), lambda x: Polygon(x.exterior).buffer(0)):
        try:
            v = make_valid(f(g)); parts = sorted([q for q in getattr(v, "geoms", [v]) if q.geom_type == "Polygon"], key=lambda q: -q.area)
            if parts and parts[0].is_valid: return parts[0]
        except Exception: pass
    return g
