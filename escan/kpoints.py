""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ['KPoints', 'KGrid', 'ReducedKGrid', 'ReducedKDensity']
from abc import ABCMeta, abstractmethod
 
class KPoints(object):
  """ Abstract base class for callable KMesh objects. 

      It is expected that the kpoints are expressed in cartesian coordinates
      of the reciprocal space. There is no 2|pi|. In other words, the
      following code holds true:

      >>> for count, kvec in kpoints_instance(structure): 
      >>>   assert abs(dot(kvec, structure.cell[:,0]) * 2e0 * pi -1e0) < 1e-12

      The code above assumes that structure is a valid `crystal.Structure`
      object.

      .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
  """
  __metaclass__ = ABCMeta

  def __init__(self): object.__init__(self)

  @abstractmethod
  def _mnk(self, input, output):
    """ Returns iterable yielding (multiplicity, kpoints). 
    
        :Parameters:
          input : `lada.crystal.Structure`
            The structure before vff relaxation (if any). 
          output : `lada.crystal.Structure`
            The structure after vff relaxation. If no relaxation was performed,
            then should be the same structure as the input structure.
    """
    pass
  @abstractmethod
  def __repr__(self):
    """ This object must be representable. """
    pass
  
  def multiplicity(self, input, output):
    """ Generator yielding multiplicity. """
    for count, k in self._mnk(input, output): yield count
  def kpoint(self, input, output):
    """ Generator yielding kpoints. """
    for count, k in self._mnk(input, output): yield k
  def __call__(self, input, output): 
    """ Iterator over (multiplicity, kpoint) tuples. """
    for r in self._mnk(input, output): yield r
    

class KContainer(object):
  """ Simple KPoints class which acts as a container. """
  def __init__(self, kpoints, multiplicity, relax=True):
    """ Initializes the kpoint container. """
    self.kpoints = [k for k in kpoints]
    """ Sequence of kpoints. """
    self.multiplicity = multiplicity
    """ Sequence with the multiplicity of the respective kpoints. """
    if self.multiplicity == None:
      self.multiplicity = [1e0 / len(self.kpoints) for k in self.kpoints]
    else: self.multiplicity = [m for m in self.multiplicity]
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  def _mnk(self, input, output):
    """ Loop over array of kpoints. 
    
        If ``self.relax`` is True, then kpoints are relaxed from the input cell
        to the relaxed output cell.
    """
    from numpy import dot, abs
    from numpy.linalg import inv

    if self.relax: deformation = dot(inv(output.cell.T), input.cell.T)
    for count, k in zip(self.kpoints, self.multiplicity):
      if self.relax: k = dot(deformation, k)
      yield count, k
    
  def kpoint(self, input, output):
    """ Generator yielding kpoints. """
    for count, k in self._mnk(input, output): yield k
  def __call__(self, input, output): 
    """ Iterator over (multiplicity, kpoint) tuples. """
    for r in self._mnk(input, output): yield r

  def __repr__(self, *args):
    return '{0.__class__.__name__}({1},{2})'\
           .format(self, repr(self.kpoints), repr(self.multiplicity))

class KGrid(KPoints):
  """ Unreduces kpoint grid with offsets. """

  def __init__(self, grid = None, offset = None, relax = True):
    """ Initializes unreduced KGrid. """
    from numpy import array
    KPoints.__init__(self)
    self.grid = grid if grid != None else array([1,1,1])
    """ Grid dimensions in reciprocal space. """
    self.offset = offset if offset != None else array([0,0,0])
    """ Offset from Gamma of the grid. """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """


  def _mnk(self, input, output):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv, norm
    from numpy import zeros, array, dot
    invcell = output.cell.T
    cell = inv(invcell)
    
    if self.relax: deformation = dot(cell, input.cell.T)

    a = zeros((3,), dtype='float64')
    weight = 1e0 / float(self.grid[0] * self.grid[1] * self.grid[2])
    for x in xrange(self.grid[0]):
      a[0] = float(x + self.offset[0]) / float(self.grid[0])
      for y in xrange(self.grid[1]):
        a[1] = float(y + self.offset[1]) / float(self.grid[1]) 
        for z in xrange(self.grid[2]):
          a[2] = float(z + self.offset[2]) / float(self.grid[2])
          yield 1, (dot(deformation, a) if self.relax else a)
          
          
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_one = all( abs(array(self.grid)-array([1,1,1])) < 1e-12 )
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_one and is_zero:
      return '{0.__class__.__name__}(relax={0.relax})'.format(self)
    if is_one:
      return '{0.__class__.__name__}(relax={0.relax}, '\
             'offset=({0.offset[0]},{0.offset[0]},{0.offset[0]}))'\
             .format(self)
    if is_zero:
      return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}), '\
             'relax={0.relax}))'.format(self)
    return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}), '\
           '({0.offset[0]},{0.offset[0]},{0.offset[0]}), relax={0.relax})'.format(self)


class KDensity(KGrid):
  """ Unreduced kpoint grid parameterized by the density and offset. """

  def __init__(self, density, offset = None, relax=True):
    """ Initializes unreduced KGrid. """
    from numpy import array
    KGrid.__init__(self, relax=relax, offset=offset)
    self.density = density
    """ 1-dimensional density in cartesian coordinates (1/Angstrom). """

  def _mnk(self, input, output):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv, norm
    from numpy import zeros, array, dot, floor
    from ..crystal import fill_structure
    from ..crystal.gruber import Reduction
    from quantities import angstrom
    
    reduction = Reduction()
    cell = reduction(output.cell, recip=True) * output.scale
    density = self.density
    if hasattr(density, 'rescale'): density.rescale(1e0/Angstrom)
    self.grid = [0,0,0]

    for i in range(3):
      self.grid[i] = norm(cell[:,i]) / self.density
      self.grid[i] = int(max(1, floor(grid[i]+0.5)))

    return KGrid._mnk(self, input, output)
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_zero:
      return '{0.__class__.__name__}({0.density}, relax={0.relax})'.format(self, repr(self.cell))
    return '{0.__class__.__name__}(({1}, ({2[0]},{2[0]},{2[0]}), relax={0.relax})'\
           .format(self, repr(self.cell), self.offset)

def _reduced_grids_factory(name, base):
  class ReducedKGrid(base): 
    def __init__(self, *args, **kwargs):
      """ Initializes reduces k-grid.
      
          :param args: Passed on to base class.
          :param kwargs: Passed on to base class.
          :kwarg tolerance: real number which defines the criterion by which
            two vectors are recognized as equal. Defaults to 1e-12.
      """
      self.tolerance = kwargs.pop('tolerance', 1e-12)
      """ Criteria to determine whether two k-vectors are the same. """
      base.__init__(self, *args, **kwargs)

    def _mnk(self, input, output):
      """ Returns list of inequivalent vectors with multiplicity. """
      from numpy.linalg import inv, norm
      from ..crystal import Lattice, to_origin, to_voronoi
      recip = Lattice()
      recip.cell = inv(output.cell.T)
      recip.add_site = (0,0,0), '0'
      recip.find_space_group()

      # now checks whether symmetry kpoint exists or not.
      seen = []
      for mult, kpoint in base._mnk(self, input, output): 
        found = False
        kpoint = to_origin(kpoint, recip.cell)
        for i, (count, vec) in enumerate(seen):
          for op in recip.space_group:
            u = to_voronoi(op(kpoint), recip.cell, vec)
            if all(abs(u) < self.tolerance):
              found = True
              seen[i][0] += mult
              break
          if found: break
        if found == False: seen.append([mult, kpoint.copy()])
      return seen

    def unreduced(self, input, output):
      """ Yields unreduced kpoints. """
      for mult, kpoint in base._mnk(self, input, output): yield mult, kpoint

    def mapping(self, input, output):
      """ Yields index of unreduced kpoint in array of reduced kpoints. """
      from numpy.linalg import inv, norm
      from ..crystal import Lattice, to_origin, to_voronoi
      recip = Lattice()
      recip.cell = inv(output.cell.T)
      recip.add_site = (0,0,0), '0'
      recip.find_space_group()
      
      # now checks whether symmetry kpoint exists or not.
      seen = []
      for mult, kpoint in base._mnk(self, input, output): 
        found = False
        kpoint = to_origin(kpoint, recip.cell)
        for i, (count, vec) in enumerate(seen):
          for op in recip.space_group:
            u = to_voronoi(op(kpoint), recip.cell, vec)
            if all(abs(u) < self.tolerance):
              found = True
              seen[i][0] += mult
              yield i
              break
          if found: break
        if found == False:
          seen.append([mult, kpoint.copy()])
          yield len(seen)-1

    def __repr__(self):
      """ Represents this object. """
      if self.tolerance == 1e-12: return base.__repr__(self)
      result = base.__repr__(self)
      result = result[:-1].rstrip()
      if result[-1] == '(': return result + 'tolerance={0})'.format(self.tolerance)
      return result + ', tolerance={0})'.format(self.tolerance)
  ReducedKGrid.__name__ = name
  ReducedKGrid.__doc__ = """ {0} reduced according to symmetries. """.format(base.__name__)
  ReducedKGrid.__module__ = base.__module__ 
  return ReducedKGrid


ReducedKGrid    = _reduced_grids_factory('ReducedKGrid', KGrid)
ReducedKDensity = _reduced_grids_factory('ReducedKDensity', KDensity)
