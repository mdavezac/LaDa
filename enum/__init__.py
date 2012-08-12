__docformat__ = "restructuredtext en"

class Transforms(object):
  """ Lattice transformation object.

      This object can create the correct permutation for any transformation of
      the periodic lattice (except for pure translation).
  """
  def __init__(self, lattice):
    """ Creates Transform object.

        :param op:
          4x3 matrix representing the transformation.
        :param lattice: 
          Lattice which forms the back-bone on which to enumerate structures.
        :type lattice: :py:class:`~lada.crystal.cppwrappers.Structure`
    """
    from numpy import dot
    from numpy.linalg import inv
    from ..crystal import which_site, Structure, into_voronoi, space_group
    super(Transform, self).__init__()
    self.lattice = lattice.copy()
    self.space_group = space_group(self.lattice)
    sites = []
    for site in lattice:
      if not hasattr(site.type, '__iter__'): continue
      if len(site.type) > 1: sites.append(site)
    self.dnt = []
    """ Site permutations and translation vector. """
    invcell = inv(lattice.cell)
    for op in self.lattice.space_group:
      self.dnt.append([])
      for i, site in enumerate(sites):
        newpos = dot(op[:3], site.pos) + op[3]
        j = which_sites(newpos, sites, invcell)
        trans = into_voronoi(site.pos - newpos, lattice.cell, invcell)
        self.dnt[-1].append( (j-i, trans) )

  def transformations(self, hft):
    """ Creates permutations for given Hart-Forcade transform. """
    from itertools import product
    from numpy import zeros, dot
    from numpy.linalg import inv
    result = zeros( (len(self.lattice.space_group), hft.size * len(self.tnd[0])),
                    dtype='int')
    for n, op in enumerate(self.lattice.space_group):
      rotation = dot(hft.transform, dot(op[:3], inv(hft.transform)))
      iterpos = [ xrange(hft.quotient[0]), 
                  xrange(hft.quotient[1]), 
                  xrange(hft.quotient[2]) ] 
      size = hft.size
      for siterperm, translation in self.tnd:
        for i,j,k in enumerate(product(*itertrans)):
          newpos = dot(rotation, [i, j, k]) + translation
          k = round(newpos[0]+1e-6)
          l = round(newpos[1]+1e-6)
          m = round(newpos[2]+1e-6)
          result[n, siterpem*site+self.lattice_index(i, j, k, hft)]            \
              = self.lattice_index(k, l, m, hft)
    return result

  def translations(self, hft):
    """ Array of permutations arising from pure translations """
    from itertools import product
    from numpy import zeros
    nsites = len(self.tnd[0])
    itertrans = [ xrange(hft.quotient[0]), 
                  xrange(hft.quotient[1]), 
                  xrange(hft.quotient[2]) ] 
    size = hft.size
    result = zeros((size, nsites * size), dtype='int64') 
    for n, (i,j,k) in enumerate(product(*itertrans)):
      iterpos = [ xrange(hft.quotient[0]), 
                  xrange(hft.quotient[1]), 
                  xrange(hft.quotient[2]) ] 
      for l, m, n in product(*iterpos):
        i = (i+l) % hft.quotient(0)
        j = (j+m) % hft.quotient(1)
        k = (k+n) % hft.quotient(2)
        result[n, self.lattice_index(l, m, n, hft)::nsites]                    \
            = self.lattice_index(i, j, k, hft)
    return result
  
  @staticmethod
  def lattice_index(i, j, k, hft):
    return k + (hft.quotient[2] * (j + i * hft.quotient[1])) 

