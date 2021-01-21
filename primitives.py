import math

def sqdist(p, q):
  """Squared Euclidean distance between p and q"""
  
  return (p[0]-q[0]) * (p[0]-q[0]) + (p[1]-q[1]) * (p[1]-q[1])

def dist(p, q):
  """Euclidean distance between p and q"""

  return math.sqrt(sqdist(p,q))

def signed_area(p, q, r):
  """Signed area of the triangle p,q,r"""

  return (((q[0]-p[0]) * (r[1]-p[1])) - ((r[0]-p[0]) * (q[1]-p[1]))) / 2

def area(p, q, r):
  """Area of the triangle p,q,r"""

  return abs(signed_area(p,q,r))

def ccw(p, q, r):
  """True if p,q,r are oriented counterclockwise"""

  return signed_area(p, q, r) > 0

def colinear(p, q, r):
  """True if p,q,r are colinear"""

  return signed_area(p, q, r) == 0

def bbox(v):
  """Bounding box of a set of points"""

  def get_y(p):
    return p[1]
  xmin = min(v)[0]
  ymin = min(v, key = get_y)[1]
  xmax = max(v)[0]
  ymax = max(v, key = get_y)[1]
  return (xmin,ymin), (xmax,ymax)

def point_in_rect(p, rect):
  """True if point p is inside the bounding box of rect"""

  (xmin,ymin), (xmax,ymax) = bbox(rect)
  return xmin <= p[0] <= xmax and ymin <= p[1] <= ymax

def point_in_seg(p, s):
  """True if point p is on the segment s"""

  if not colinear(p, *s):
    return False
  return point_in_rect(p, s)

def seg_intersect(s, t):
  """True if the two segments intersect"""

  s1,s2 = s
  t1,t2 = t
  
  a1 = signed_area(s1,s2,t1)
  a2 = signed_area(s1,s2,t2)
  a3 = signed_area(t1,t2,s1)
  a4 = signed_area(t1,t2,s2)
  
  if a1 == 0 or a2 == 0 or a3 == 0 or a4 == 0: # Colinear
    return point_in_seg(s1,t) or point_in_seg(s2,t) or point_in_seg(t1,s) or point_in_seg(t2,s)
  
  return (a1 > 0) != (a2 > 0) and (a3 > 0) != (a4 > 0)

def seg_proper_intersect(s,t):
  """True if the two segments not at a common vertex"""

  return (seg_intersect(s,t)
    and s[0] != t[0] and s[1] != t[0] and s[0] != t[1] and s[1] != t[1])

def poly_area(poly):
  """Area of a polygon"""

  a = 0
  for i in range(len(poly)):
    a += (poly[i-1][0] + poly[i][0]) * (poly[i-1][1] - poly[i][1])
  return abs(a/2)

def area_change(ladd, ldel):
  """Area change when we add a set of edges and remove a set of edges"""

  a = 0
  for e in ladd:
    a += (e[0][0] + e[1][0]) * (e[0][1] - e[1][1])
  for e in ldel:
    a -= (e[0][0] + e[1][0]) * (e[0][1] - e[1][1])
  return -a/2

def score(poly):
  """Challenge score of a polygon"""

  ch = convex_hull2(poly)
  ach = poly_area(ch)
  apoly = poly_area(poly)
  return apoly / ach
  
def read_points(f):
  """Read the input file"""  
  
  indexofpoint = {}
  for line in f:
    if line and line[0] != "#":
      i,x,y = (int(s) for s in line.split())
      indexofpoint[(x,y)] = i
  return indexofpoint

def write_polygon(poly, indexofpoint, f):
  """Write the polygon to the output file"""  
  
  for p in poly:
    print(indexofpoint[p], file = f)

def _halfch2(v):
  """Graham scan to compute upper hull of sorted points"""

  hull = [v[0]]
  for p in v[1:]:
    while len(hull)>=2 and signed_area(p, hull[-1], hull[-2]) > 0:
      hull.pop()
    hull.append(p)
  return hull

def convex_hull2(v):
  """Convex hull analogue that includes vertices on edges"""

  top = _halfch2(sorted(v))
  bottom = _halfch2(sorted(v, reverse=True))
  return top[0:-1] + bottom[0:-1]
