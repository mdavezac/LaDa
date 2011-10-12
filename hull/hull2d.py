
class Hull2d(object):
  """ A 2d-convex hull. """
  dimension = 2
  """ Dimension of this convex-hull. """
  def __init__(self, vertices=None, tolerance = 1e-12):
    """ Creates a planar polygon using the vertices. """
    self.tolerance = tolerance
    """ Tolerance of fuzzy comparison operators. """
    self._basis = None
    """ Basis of the plane of the convex-hull. """
    self.vertices = vertices

  @property 
  def vertices(self):
    """ Sets of vertices. 

        This is a numpy array of vertices, arranged in clockwise order.
        The planar polygon owns a copy of the vertices.

        >>> self.vertices = vertices
        >>> assert id(self.vertices) != id(vertices)

        The first line above will assert if the vertices are not planar.
    """
    return self._vertices
  
  @vertices.setter
  def vertices(self, sample):
    """ Taken from literate programming website. """
    from numpy import array, dot, concatenate
    from numpy.linalg import norm
    # case where vertices are None.
    if sample is None: self._vertices = None; return
    sample = array(sample) # copy to  a numpy array.

    # case where there are fewer vectors than dimensions: always planar.
    assert len(sample.shape) != 1, ValueError("Empty array of points.")
    assert sample.shape[0] >= 3, ValueError("Too few points to construct convex hull.")

    # checkings that the vertices are planar.
    origin = sample[0]
    endpoint = max(sample[1:], key=lambda x: norm(x-origin))
    origin = max(sample, key=lambda x: norm(x-endpoint))
    a0 = endpoint - origin; a0 /= norm(a0)
    a1 = max(sample, key=lambda x: norm(x-origin - dot(x-origin, a0)*a0))
    a1 = a1-origin - dot(a1-origin, a0)*a0
    if norm(a1) < self.tolerance: # all colinear 
      self.vertices = array([origin, endpoint])
      self._basis = concatenate(([(0,0,0)], [a0], [origin]))
      return
    a1 /= norm(a1) # normalize basis.
    # check for coplanarity
    for vec in sample:
      centered  = vec - origin
      centered -= dot(centered, a0) * a0 + dot(centered, a1) *a1
      assert norm(centered) < self.tolerance, \
          ValueError("Vertices are not coplanar: {0}.".format(centered))
    self._basis = array([a0, a1, origin])

    # now creates convex hull
    base = Hull2d.edge(origin, endpoint)
    self._vertices = Hull2d.link(self.dome(sample, base), self.dome(sample, base[::-1]))
  
  def dome(self, sample, base):
    """ Creates part of the convex hull. """
    from numpy import repeat, argmax, dot
    h, t = base
    a0 = t - h
    a1 = dot(self._basis[0], a0) * self._basis[1] - dot(self._basis[1], a0) * self._basis[0]
    dists = dot(sample-h, a1)
    outer = repeat(sample, dists > self.tolerance, 0)
    if len(outer):
      pivot = sample[argmax(dists)]
      return Hull2d.link( self.dome(outer, Hull2d.edge(h, pivot)),\
                          self.dome(outer, Hull2d.edge(pivot, t)) )
    else: return base


  def project(self, points):
    """ Projects points on 2d-plane. """
    from numpy import dot, array
    assert self._basis is not None, RuntimeError("Convex-hull has not been constructed.")
    result = []
    for p in points:
      p = p - self._basis[-1]
      result.append( (dot(p, self._basis[0]), dot(p, self._basis[1])) )
    return array(result)

  @staticmethod
  def edge(a, b):
    """ Links two edges. """
    from numpy import concatenate
    return concatenate(([a],[b]))
  @staticmethod
  def link(a, b): 
    """ Concatenates two edges. """
    from numpy import concatenate
    return concatenate((a,b[1:]))

  def __getitem__(self, index):
    """ Returns vertex at ``index``. """
    assert "_vertices" in self.__dict__, RuntimeError("Vertices do not exist.")
    return self._vertices[index]

  def __iter__(self): 
    """ Iterates over vertices. """
    assert "_vertices" in self.__dict__, RuntimeError("Vertices do not exist.")
    return self._vertices.__iter__()

  def __len__(self):
    """ Number of vertices. """
    assert "_vertices" in self.__dict__, RuntimeError("Vertices do not exist.")
    return len(self._vertices)



def _test(n=20, d=3):
  """ Tests hull. """
  from numpy import array, concatenate
  from numpy.random import normal
  a0 = normal(size=(d))
  a1 = normal(size=(d))
  origin = normal(size=(d))
  points0 = array([origin + a0*normal() + a1*normal() for i in range(n)])
  points1 = array([origin + 2.*(0.5*a0-0.5*a1)*normal() for i in range(min(n,50))])
  points = concatenate((points0, points1))
  hull = Hull2d(points)
  return hull.project(points), hull.project(hull.vertices)
  

