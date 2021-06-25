#!/usr/bin/pypy3
import sys
import math
import random
import time
import heapq
from primitives import *
from geometer import *

starttime = time.time()

# Random seed
seed = 1

# Maximize or minimize
maximize = True

# Parameter of the weight function
# Referred to as 1/alpha in the paper
pen = 90

# Run the algorithm multiple times and choose best solution
multirun = False

# Maximum time in seconds for a new execution to be started
timeout = 150

# Gaussian noise standard deviation applied to the weight function
sigma = 0

# Apply local search optimization
opt = True

# Stop if optimization gain is less than optgain
optgain = .001

# Maximum length of a path that can be moved in the local search
# Referred to as ell in the paper
hops = 1

# How far in terms of grid cells around an edge we should check for candidate points
# Number of cells is (2*hood + 1)**2
# Referred to as kappa in the paper
hood = 2

# Abort if the number of points is smaller than nmin
nmin = 0

# Abort if the number of points is greater than nmax
nmax = 100000

comments = ""

def perturb(x):
  """Perturb sigma with Gaussian noise"""
  
  return x + x * abs(random.gauss(0.0, sigma))

def greedy_start(pts):
  """Run the greedy algorithm from a set of points"""
  if not multirun: print("Started greedy solver")
  
  if maximize:
    # Start with the convex hull including vertices on edges
    poly = convex_hull2(pts)
  else:
    # Start with a random triangle
    p1 = random.choice(pts)
    p2 = min(set(pts) - {p1}, key = lambda p : sqdist(p,p1))
    p3 = min(set(pts) - {p1,p2}, key = lambda p : dist(p,p1) + dist(p,p2))
    # Make sure area is negative
    if ccw(p1,p2,p3):
      poly = [p2,p1,p3]
    else:
      poly = [p1,p2,p3]
    
  geo = Geometer(pts)
  geo.add_poly(poly)

  return greedy(geo)
  
def weight(p1,p2,p3):
  """Weight of triangle p1p2p3"""
  
  a = signed_area(p1,p2,p3)
  a += (sqdist(p1,p2) - sqdist(p2,p3) + sqdist(p3,p1)) / pen
  if sigma:
    a = perturb(a)
  return a

def greedy(geo, nhood = None):
  """Run the greedy algorithm from a geometer object"""
  if not nhood:
    nhood = hood
  
  todo = geo.pts - set(geo.get_poly())
  dic_ep = {}
  cand = []
  
  def create_dic(e):
    if nhood == math.inf:
      cand_pts = todo
    else:
      cand_pts = todo & geo.pts_near_seg(e, nhood)
    
    if maximize:
      dic_ep[e] = [(weight(p,*e),p,e) for p in cand_pts]
    else:
      dic_ep[e] = [(weight(p,*e),p,e) for p in cand_pts
                                      if signed_area(p,*e) > 0]
    dic_ep[e].sort(reverse = True)
  
  def get_point(e):
    while dic_ep[e]:
      ape = dic_ep[e][-1]
      if ape[1] in todo:
        return ape
      dic_ep[e].pop()
    return (math.inf,None,e)
      
  for e in geo:
    create_dic(e)

  while todo:
    if not multirun and len(todo) % int(1+n/50) == 0:
      print(".", end="")
      sys.stdout.flush()
    
    cand.clear()
    for e in dic_ep:
      cand.append(get_point(e))
    heapq.heapify(cand)

    while cand:
      a,p,e = cand[0]
      if not p:
        if nhood != math.inf:
          # No solution will be found! Try to switch to infinite neighborhood
          if not multirun: print("\nSetting hood to inf with",len(todo),"points left")
          return greedy(geo, math.inf)
        else:
          return None
      add_list = [(e[0],p),(p,e[1])]
      del_list = [e]
      if geo.valid_change(add_list, del_list):
        # Found valid triangle
        todo.remove(p)
        geo.apply_change(add_list, del_list)
        del dic_ep[e]
        for ne in add_list:
          create_dic(ne)
        break
      dic_ep[e].pop()
      heapq.heappushpop(cand, get_point(e))

  if not multirun: print("")

  return geo.get_poly()

def local_search_step(poly):
  """Returns the local search optimized polygon"""
  
  geo = Geometer(poly)
  geo.add_poly(poly)
  n = len(poly)
  
  def gen_epath():
    for e in geo:
      if hood == math.inf:
        vertices = geo.vertices()
      else:
        vertices = geo.pts_near_seg(e,hood)
      for p1 in vertices:
        path = []
        p = p1
        for h in range(hops):
          if p in e:
            break
          path.append(p)
          yield (e,path)
          p = geo.next_vertex[p]
  
  changes_todo = []
 
  for e,path in gen_epath():
    before_path = geo.prev_vertex[path[0]]
    after_path = geo.next_vertex[path[-1]]
    add_list = [(e[0],path[-1]), (path[0],e[1]), (before_path,after_path)]
    del_list = [(before_path,path[0]), (path[-1],after_path), e]
    for i in range(1,len(path)):
      add_list.append((path[i],path[i-1]))
      del_list.append((path[i-1],path[i]))

    a = area_change(add_list, del_list)
    if a > 0 and geo.valid_change(add_list, del_list):
      changes_todo.append((a, add_list, del_list))

  changes_todo.sort(reverse=True)
  
  for _,add_list,del_list in changes_todo:
    if geo.valid_change(add_list, del_list):
      geo.apply_change(add_list, del_list)

  return geo.get_poly()

def local_search(poly):
  """Runs local search until score does not change much"""
  if not poly:
    return None

  if not multirun: print("Started refining with:", score(poly))

  delta = None

  while delta == None or delta >= optgain:
    oldscore = score(poly)
    poly = local_search_step(poly)
    delta = abs(oldscore - score(poly))
    if not multirun: print("Improved to", score(poly))
  return poly

def manyruns(pts):
  """Runs greedy + local search until timeout"""
  endtime = starttime + timeout
  bestsol = None
  global comments
  
  if not multirun:
    poly0 = greedy_start(pts)
    if opt:
      comments += "# Time before opt: " + str(time.time()-starttime) + "\n"
      comments += "# Score before opt: " + str(score(poly0)) + "\n"
      print(comments)
      bestsol = local_search(poly0)
    else:
      bestsol = poly0

  else:
    while True:
      poly0 = greedy_start(pts)
      poly = local_search(poly0)
      
      if poly and (not bestsol or
                   (maximize and poly_area(poly) > poly_area(bestsol)) or
                   (not maximize and poly_area(poly) < poly_area(bestsol))):
        print(f"{int(time.time() - starttime)} sec, {score(poly0)} => {score(poly)}")
        bestsol = poly

      if time.time() > endtime:
        break
    
  return bestsol


def save(sol):
  """Saves the solution sol to a file.
  The file name is determined based on global variables"""
  
  extension = ""
  if maximize:
    extension += "max"
  else:
    extension += "min"
  if pen != 90:
    extension += "pen" + str(pen)
  if sigma:
    extension += "sigma" + str(int(100*sigma))
  if hood != math.inf:
    extension += "hood" + str(hood)
  if opt:
    extension += "opt"
    if hops > 1:
      extension += str(int(hops))

  fn = basename+"." + extension + ".solution"
  print("Writing", fn)

  with open(fn, "w") as outfile:
    print("# Score:", score(sol), file = outfile)
    print("# Time:", time.time()-starttime, file = outfile)
    print("# Parameters:",sys.argv, file = outfile)
    print(comments, file = outfile,end="")  
    write_polygon(sol, indexofpoint, outfile)


#---------------------------------------------------------
# Actual execution starts here
#---------------------------------------------------------

if len(sys.argv) < 2:
  print("""
pypy3 poLYG.py [args] inputfile

Inputfile extensions are ignored. Args are of the form variable=value (no space, no dash). In fact, they will be executed with python exec(arg), so feel free to add formulas that depend on n or other parameters.  Some variables and their default values are:
""")
  print(f"maximize={maximize}   % Maximum or minimum area")
  print(f"pen={pen}   % Parameter 1/alpha of the weight function")
  print(f"hood={hood}   % Neighborhood kappa of an edge") 
  print(f"opt={opt}   % Apply local search optimization") 
  print(f"hops={hops}   % Value of ell for the local search") 
  print(f"multirun={multirun}   % Run many times")
  print(f"sigma={sigma}   % Gaussian noise of the weight function")
  print(f"seed={seed}   % Random seed")
  print(f"timeout={timeout}   % Maximum number of seconds to start a new run")
  exit()

basename = sys.argv[-1].split(".")[0]

with open(basename+".instance") as inputfile:
  indexofpoint = read_points(inputfile)

n = len(indexofpoint)

for arg in sys.argv[1:-1]:
  exec(arg)

if n > nmax:
  print("File is too large:",n)
  exit()

if n < nmin:
  print("File is too small:",n)
  exit()

random.seed(seed)
print("\n----------",basename,"started")
  
sol = manyruns(list(indexofpoint.keys()))

if not sol:
  print("No solution found!")
  exit()

print("After",time.time()-starttime, "sec on", basename, "got", score(sol))
save(sol)

