import sys
from primitives import *
import itertools

# Cell of segment s

class Geometer:
  
  def __init__(self, pts):
    self.grid_edges = {}
    self.next_vertex = {}
    self.prev_vertex = {}
    self.long_edges = set()
    self.doublearea = 0.0

    self.pts = set(pts)
    self.box = bbox(pts)
    (xmin,ymin),(xmax,ymax) = self.box
    ncols = 4 * len(pts)**.25
    self.cell_size = max(xmax-xmin,ymax-ymin) / ncols
    self.cell_size = 2 + 2 * int(self.cell_size / 2)
    self.grid_pts = {}
    for p in self.pts:
      c = self.cell(p)
      if c in self.grid_pts:
        self.grid_pts[c].add(p)
      else:
        self.grid_pts[c] = {p}

  def cell(self, p):
    return (int(p[0] / self.cell_size), int(p[1] / self.cell_size))

  def cell_box(self, c):
    return ((c[0] * self.cell_size, c[1] * self.cell_size), ((c[0]+1) * self.cell_size, (c[1]+1) * self.cell_size))
  
  def cells_near_cell(self, c, delta):
    ret = set()
    for i in range(-delta,delta+1):
      for j in range(-delta,delta+1):
        ret.add((c[0]+i, c[1]+j))
    return ret

  def cells_near_seg(self, seg, delta):
    ret = set()
    for c in self.cells(seg):
      ret |= self.cells_near_cell(c, delta)
    return ret
    
  def pts_near_seg(self, seg, delta):
    if len(self.grid_pts) <= (1+2*delta)**2:
      return self.pts
    
    ret = set()
    for c in self.cells_near_seg(seg, delta):
      if c in self.grid_pts:
        ret |= self.grid_pts[c]
    return ret
    
  def cells(self, e):
    if e[0] > e[1]:
      e = (e[1],e[0])

    c0 = self.cell(e[0])
    c1 = self.cell(e[1])
    if  c0 == c1:
      return [c0]

    upwards = e[1][1] > e[0][1]

    l = [c0]
    while l[-1] != c1:
      c = l[-1]
      (xmin,ymin), (xmax,ymax) = self.cell_box(c)
      if   seg_intersect(e, ((xmax,ymin), (xmax,ymax))):
        c = (c[0]+1, c[1])
      elif upwards:
        c = (c[0], c[1]+1)
      else:
        c = (c[0], c[1]-1)
        
      l.append(c)
      
    return l

  def add(self, e):
    self.next_vertex[e[0]] = e[1]
    self.prev_vertex[e[1]] = e[0]
    
    self.doublearea += (e[0][0] + e[1][0]) * (e[0][1] - e[1][1])

    l = self.cells(e)
    if len(l) <= 4:
      for c in l:
        if c in self.grid_edges:
          self.grid_edges[c].add(e)
        else:
          self.grid_edges[c] = {e}
    else:
      self.long_edges.add(e)
    
  def remove(self, e):
    del self.next_vertex[e[0]]
    del self.prev_vertex[e[1]]

    self.doublearea -= (e[0][0] + e[1][0]) * (e[0][1] - e[1][1])
    
    l = self.cells(e)
    if len(l) <= 4:
      for c in l:
        self.grid_edges[c].remove(e)
        if not self.grid_edges[c]:
          del self.grid_edges[c]
    else:
      self.long_edges.remove(e)

  def add_poly(self, poly):
    for i in range(len(poly)):
      e = (poly[i-1], poly[i])
      self.add(e)
      
    if self.long_edges:
      # Sets in pypy3 are sorted. Keep longer edges first
      l = list(self.long_edges)
      l.sort(key = lambda e : sqdist(*e), reverse = True)
      self.long_edges = set(l)

  def get_poly(self):
    poly = [next(iter(self.next_vertex.keys()))]
    while True:
      p = self.next_vertex[poly[-1]]
      if p == poly[0]:
        break
      poly.append(p)

    return poly
    
  def intersections(self, e):
    for s in self.long_edges:
      if seg_intersect(s, e):
        yield s
      
    for c in self.cells(e):
      for s in self.grid_edges.get(c,[]):
        if seg_intersect(s, e):
          yield s

  def valid_change(self, ladd, ldel):
    # Check if edges to delete are present
    for e in ldel:
      if e not in self:
        return False

    # Check if there are no duplicates
    if len(set(ldel)) != len(ldel) or len(set(ladd)) != len(ladd):
      return False

    # Check if there are no loop edges
    for e in ladd:
      if e[0]==e[1]:
        return False

    # Check if there are no mutual intersections
    for e1,e2 in itertools.combinations(ladd,2):
        if seg_proper_intersect(e1,e2):
          return False

    # Check if area does not change sign (important for minimization)
    newarea = self.doublearea + area_change(ladd,ldel)
    if (newarea > 0 and self.doublearea < 0) or (newarea < 0 and self.doublearea > 0):
      return False
        
    # Check if there are no intersections with polygon
    for e in ladd:
      if (e[1],e[0]) not in self and self.proper_intersects(e):
        return False
    return True

  def apply_change(self, ladd, ldel):
    for e in ldel:
      self.remove(e)
    for e in ladd:
      self.add(e)

  def __contains__(self, e):
    return self.next_vertex.get(e[0], None) == e[1]

  def __iter__(self):
    return iter(self.next_vertex.items())
      
  def vertices(self):
    return self.next_vertex.keys()

  def intersects(self, e):
    return any(self.intersections(e))

  def proper_intersects(self, e):
    for seg in self.intersections(e):
      if seg[0] != e[0] and seg[1] != e[0] and seg[0] != e[1] and seg[1] != e[1]:
        return True
    return False
  
