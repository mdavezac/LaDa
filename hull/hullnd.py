class Facet(object):
  """ Represents a facet on a hyperplane. """
  def __init__(self, vertices, tolerance=1e-12):
    """ Initializes facet. 
    
        This cannot be anything but a true hyperplane: the vertices cannot be
        linearly correlated.
    """
    self.tolerance = tolerance
    """ Tolerance of fuzzy comparison operators. """
    self.dimension = len(vertices) - 1
    """ Dimension of the triangle """
    self.vertices = vertices
    """ Vertices of the facet. """

  @property
  def vertices(self): 
    """ All points leading to this facet. """
    return self._vertices 
  @vertices.setter
  def vertices(self, sample):
    """ Construct vertices and edges from sample.
    
    
        :Parameters: 
          sample 
            An array of non-colinear points.
        
        :raise ValueError: if some points are colinear.
    """
    from itertools import combinations
    from hull3d import create_basis

    if sample is None: return
    assert len(sample.shape) != 1, ValueError("Sample too small.")
    assert sample.shape[0] == self.dimension + 1, ValueError("Sample should contain 3 vectors.")

    # check that points are non-colinear
    b = create_basis(sample, self.tolerance) 
    assert b.shape[0] == self.dimension + 1

    # sets vertices and dimension.
    self._vertices = sample

  def iteredges(self):
    """ Iterates over all edges. """
    from numpy import concatenate
    assert "_vertices" in self.__dict__, RuntimeError("Vertices not set.")
    for i in range(self.dimension+1): 
      r0 = range(i, min(self.dimension+1, self.dimension+i))
      r1 = range(0, i-1)
      vertices = self._vertices[r0 + r1]
      yield vertices

  def __eq__(self, other):
    """ Compares two facets. """
    from numpy import norm
    if self.dimension != other.dimension: return False
    for vertex in self.vertices: 
      if min([norm(u-vertex) for u in other.vertices]) > self.tolerance: return False
    return True

  def inverted(self):
    """ Inverts orientation. """
    return Facet(self._vertices[::-1], self.tolerance)

  def volume(self, vertex):
    """ Whether a vertex can be seen by this facet. """
    from numpy import ones, dot, concatenate
    from numpy.linalg import det
    d = concatenate((self._vertices, [vertex]))
    d = concatenate((d.T, [ones(shape=(d.shape[0],), dtype=d.dtype)]))
    result = det(d)
    return 0 if result < self.tolerance and result > -self.tolerance else result

class HullNd(object): 
  """ Hold/Creates a convex-hull. """
  def __init__(self, sample, dimension = 0, tolerance = 1e-12):
    """ Initializes convex-hull. """
    from numpy import array
    from hull3d import create_basis

    sample = array(sample) # copy to  a numpy array.
    self.dimension = sample.shape[1]
    """ Dimension of the space where the convex-hull is performed. """

    assert len(sample.shape) != 1, ValueError("Empty array of points.")
    assert sample.shape[0] > self.dimension, ValueError("Too few points.")

    self.tolerance = tolerance
    """ Tolerance of fuzzy comparison operators. """

    self.basis = create_basis(sample, self.tolerance)
    assert self.basis.shape[0] == self.dimension + 1, ValueError("Not a 3d space")

    self._create_facets(sample)
    """ Vertices comprising the convex-hull. """

  @property
  def facets(self):
    """ Sets of facets and vertices. """
    return self._facets
  
  def _create_facets(self, sample):
    """ Creates convex-hull from sample. """
    if sample is None: return;
    base = self._start(sample)
    self._facets = [base]
    self._dome(sample, self._facets[-1])
    self._facets.append(base.inverted())
    self._dome(sample, self._facets[-1])

  def _start(self, sample):
    """ Computes starting base. 
    
        Base consists of first d vertices which are on the hull. Vertices which
        are minimum or maximum for any one direction must be on the hull.
    """
    from numpy import argmin, argmax, ones, dot
    from hull3d import Triangle
    origin = self.basis[-1]

    # finds index of convex-hull points first.
    base = []
    for a in self.basis[:-1]:
      dists = [dot(x-origin, a) for x in sample]
      # finds minimum
      mini = argmin(dists)
      if mini not in base: base.append(mini)
      # finds maximum.
      maxi = argmax(dists)
      if maxi not in base: base.append(maxi)

    Type = Triangle if self.dimension == 3 else Facet
    assert len(base) >= self.dimension, RuntimeError("Could not find starting base.")
    return Type(sample[base[:self.dimension]], self.tolerance)

  def _dome(self, sample, base):
    """ Creates one side of the convex hull. """
    from numpy import ones, dot, repeat, argmax, concatenate, array
    from numpy.linalg import det
    from hull3d import Triangle 

    # now creates list of volumes.
    dists = array([base.volume(u) for u in sample])

    # keep only positive volumes for next iteration.
    outer = sample[dists > self.tolerance]
    # go to next step.
    if len(outer):
      pivot = sample[argmax(dists)]
      volumes = [facet.volume(pivot) for facet in self._facets]
      deleting = [facet for facet, vol in zip(self._facets, volumes) if vol >= -self.tolerance]
      self._facets = [facet for facet, vol in zip(self._facets, volumes) if vol < -self.tolerance]
      Type = Triangle if self.dimension == 3 else Facet
      for edge in self._horizon(deleting):
        self._facets.append( Type(concatenate((edge, [pivot]))) )
      for facet in self._facets: self._dome(outer, facet)

  def _horizon(self, deleting):
    """ Iterates over horizon edges. """
    from numpy import min, all
    def compare(a, b):
      for v0 in a: 
        if all( [any(v0-v1) > self.tolerance for v1 in b] ): return False
      return True
    result = []
    for facet in deleting:
      for edge in facet.iteredges():
        found = False
        for i, edge2 in enumerate(result):
          if compare(edge2, edge): found = True; break
        if found: result.pop(i)
        else: result.append(edge)
    return result

  def vertices(self):
    """ Vertices on the convex-hull. """
    result = []
    for facet in self.facets:
      for v in facet.vertices:
        if all([any(abs(u-v)) > self.tolerance for u in result]): 
          result.append(v)
    return result

  def mesh(self):
    """ Returns triangular mesh data for mayavi. """
    from numpy import array
    vertices = [(i*3, i*3+1, i*3+2) for i, f in enumerate(self.facets)]
    result = [u for facet in self.facets for u in facet.vertices]
    result = [u for u in array(result).T]
    result.append(vertices)
    return result
 

def _test(n=5, d=3, tolerance=1e-12):
  """ Tests hull. """
  from numpy import array, concatenate, zeros, diag, ones
  from numpy.linalg import norm
  from numpy.random import normal
  from hull2d import Hull2d
  points = array([normal(size=d) for i in range(n)])
  hull = HullNd(points)
  for facet in hull.facets:
    print all(facet.volume(u) <= tolerance for u in points)
  hullvertices = hull.vertices()
  for n in range(d):
    indices = range(d)
    indices.pop(n)
    onhull = [u[indices] for u in hullvertices]
    if d == 3: 
      nm1 = Hull2d([u[indices] for u in points])
      nm1vertices = nm1.vertices
    else: 
      nm1 = HullNd([u[indices] for u in points])
      nm1vertices = nm1.vertices()
    for v0 in nm1vertices:
      found = False
      assert any(all(abs(v1-v0) < tolerance) for v1 in onhull)
  return hull, points
