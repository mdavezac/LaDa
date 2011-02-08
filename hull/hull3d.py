def create_basis(sample, tolerance=1e-12):
  """ Creates an othornormal basis.
  
      :Parameters: 
        sample 
          Array of points. Should be all of the same dimension. There should be
          more than one point
        tolerance 
          Fuzzy comparison.

      :return: A numpy array where the last (row) vector is the origin, and the
        others form an orthonormal basis for sample.
  """
  from numpy import dot, array
  from numpy.linalg import norm, det

  sample = array(sample)
  assert len(sample.shape) > 1, RuntimeEror("Cannot construct basis from a single point.")

  origin = sample[0]
  endpoint = max(sample[1:], key=lambda x: norm(x-origin))
  origin = max(sample, key=lambda x: norm(x-endpoint))
  basis = [(endpoint-origin) / norm(endpoint-origin)]

  def perpendicular(current):
    """ Returns vector perpendicular to basis. """
    result = current - origin
    for a in basis: result -= dot(result, a) * a
    return result

  while True:
    perps = [perpendicular(v) for v in sample]
    next = max(perps, key = lambda x: norm(x))
    if norm(next) < tolerance: break
    basis.append(next / norm(next))

  if len(basis) == basis[0].shape[0] and det(basis) < 0 and len(basis) > 2:
    basis[-1], basis[-2] = basis[-2], basis[-1]

  basis.append(origin)
  return array(basis)
    

class Triangle(object):
  """ Represents a facet on a hyperplane. """
  dimension = 2
  def __init__(self, vertices, tolerance=1e-12):
    """ Initializes facet. 
    
        This cannot be anything but a true hyperplane: the vertices cannot be
        linearly correlated.
    """
    self.tolerance = tolerance
    """ Tolerance of fuzzy comparison operators. """
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

    if sample == None: return
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
    yield concatenate(([self._vertices[0]], [self._vertices[1]]))
    yield concatenate(([self._vertices[1]], [self._vertices[2]]))
    yield concatenate(([self._vertices[2]], [self._vertices[0]]))

  def __eq__(self, other):
    """ Compares two facets. """
    from numpy import norm
    if self.dimension != other.dimension: return False
    for vertex in self.vertices: 
      if min([norm(u-vertex) for u in other.vertices]) > self.tolerance: return False
    return True

  def inverted(self):
    """ Inverts orientation. """
    return Triangle(self._vertices[::-1], self.tolerance)

  def _basis(self):
    """ Creates a basis from triangle. """
    from numpy import concatenate, dot, cross
    from numpy.linalg import norm
    a0 = self._vertices[1] - self._vertices[0]; a0 /= norm(a0)
    a1 = self._vertices[2] - self._vertices[1]; a1 -= dot(a1, a0); a1 /= norm(a0)
    a2 = cross(a0, a1)
    return concatenate(([a0], [a1], [a2], [self._vertices[0]]))


  def volume(self, vertex):
    """ Whether a vertex can be seen by this facet. """
    from numpy import ones, dot, concatenate
    from numpy.linalg import det
    basis = self._basis()
    d = ones((self.dimension+2, self.dimension+2), dtype="float64")
    for i, a in enumerate(basis[:-1]):
      for j, v in enumerate(self.vertices): 
        d[i,j] = dot(v-basis[-1], a) 
      d[i,-1] = dot(vertex-basis[-1], a)
    result = det(d)
    return 0 if result > -self.tolerance and result < self.tolerance else result

  def __str__(self):
    """ Prints triangle. """
    return "Triangle: {0._vertices[0]}, {0._vertices[1]}, {0._vertices[2]}.".format(self)

class Hull3d(object): 
  """ Hold/Creates a convex-hull. """
  dimension = 3
  def __init__(self, sample, dimension = 0, tolerance = 1e-12):
    """ Initializes convex-hull. """
    from numpy import array
    sample = array(sample) # copy to  a numpy array.
    assert len(sample.shape) != 1, ValueError("Empty array of points.")
    assert sample.shape[0] > self.dimension, ValueError("Too few points.")
    # initializes values.
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
    if sample == None: return;
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

    assert len(base) >= self.dimension, RuntimeError("Could not find starting base.")
    return Triangle(sample[base[:self.dimension]], self.tolerance)

  def _dome(self, sample, base):
    """ Creates one side of the convex hull. """
    from numpy import ones, dot, repeat, argmax, concatenate, array
    from numpy.linalg import det

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
      for edge in self._horizon(deleting):
        self._facets.append( Triangle(concatenate((edge, [pivot]))) )
      for facet in self._facets: self._dome(outer, facet)

  def _horizon(self, deleting):
    """ Iterates over horizon edges. """
    from numpy import min, all
    def compare(a, b):
      if all( abs(a[0]-b[0]) < self.tolerance ):
        return all( abs(a[1]-b[1]) < self.tolerance )
      if all( abs(a[0]-b[1]) < self.tolerance ):
        return all( abs(a[1]-b[0]) < self.tolerance )
      return False
    result = []
    for facet in deleting:
      for edge in facet.iteredges():
        found = False
        for i, edge2 in enumerate(result):
          if compare(edge2, edge): found = True; break
        if found: result.pop(i)
        else: result.append(edge)
    return result

  def mesh(self):
    """ Returns triangular mesh data for mayavi. """
    from numpy import array
    x, y, z, vertices = [], [], [], []
    for i, facet in enumerate(self.facets):
      x.extend(facet.vertices[:, 0])
      y.extend(facet.vertices[:, 1])
      z.extend(facet.vertices[:, 2])
      vertices.append((i*3, i*3+1, i*3+2))
    return x, y, z, vertices
 

def _test(n=5):
  """ Tests hull. """
  from numpy import array, concatenate, zeros
  from numpy.linalg import norm
  from numpy.random import normal
  d = 3
  a0 = zeros((3,), dtype="float64"); a0[0] = 1
  a1 = zeros((3,), dtype="float64"); a1[1] = 1
  a2 = zeros((3,), dtype="float64"); a2[2] = 1
  points = array([a0*normal() + a1*normal() + a2*normal() for i in range(n)])
  hull = Hull3d(points)
  return hull, points
  

