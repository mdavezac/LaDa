""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ["KEscan", 'KPoints', 'KGrid', 'ReducedKGrid', 'ReducedKDensity']
from abc import ABCMeta, abstractmethod
from .functional import Functional
from .. import __all__ as all_lada_packages
 
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
  def _mnk(self, structure):
    """ Returns iterable yielding (multiplicity, kpoints). """
    pass
  @abstractmethod
  def __repr__(self, structure):
    """ This object must be representable. """
  
  def multiplicity(self, structure):
    """ Generator yielding multiplicity. """
    for count, k in self._mnk(structure): yield count
  def kpoint(self, structure):
    """ Generator yielding kpoints. """
    for count, k in self._mnk(structure): yield k
  def __call__(self, structure): 
    """ Iterator over (multiplicity, kpoint) tuples. """
    for r in self._mnk(structure): yield r
    

class KEscan(Functional):
  """ A wrapper around Escan for computing many kpoints. """
  def __init__(self, kpoints = None, multiplicity=None):
    """ Initializes the KEscan functional. """
    self.kpoints = kpoints
    """ Kpoints to use for calculations.
    
        This object must be None (Gamma), a KMesh-derived structure, a single
        kpoint, or a list of kpoints.

        It is expected that the kpoints are expressed in cartesian coordinates
        of the reciprocal space. There is no 2|pi|. See `Kpoints` for more
        information.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """

  # need jobs package to run this code.
  if 'jobs' in all_lada_packages: 
    def __call__(self, structure, outdir=None, comm=None, **kwargs):
      """ Performs calculcations. """
      from os.path import join
      from ...jobs import JobDict, Bleeder

      kpoints = kwargs.pop('kpoints', self.kpoints)
      dont_deform_kpoints = getattr(kpoints, 'dont_deform_kpoint', self._dont_deform_kpoint)
      dont_deform_kpoints = kwargs('_dont_deform_kpoint', kpoints)
      is_mpi = False if comm == None else comm.size > 1
      is_root = True if not is_mpi else comm.rank == 0
      vffrun = kwargs.get("vffrun", None)
      genpotrun = kwargs.get("genpotrun", None)

      # performs vff calculations
      if vffrun != None: 
        vffout = super(KEscan, self)(structure, outdir, comm, do_escan=False, **kwargs)
        kwargs['vffrun'] = vffout
        if kwargs.get('genpotrun', None) == None: kwargs['genpotrun'] = vffout
  
      # create list of kpoints.
      kpoints = self._interpret_kpoints(kpoints, vffout)
      # checks for 
      if len(kpoints) == 1:
        kwargs['kpoint'] = kpoints[0]
        result = super(KEscan, self).__call__(structure, outdir, comm, *args, **kwargs)
        return result

      jobdict = JobDict()
      for i, kpoint in enumerate(kpoints):
        job = jobdict / 'kpoint_{0}'.format(kpoint)
        job.functional = self
        job.jobparams = kwargs.copy()
        job.jobparams['kpoint'] = kpoint
        job.jobparams['structure'] = structure
        job.jobparams['outdir'] = join(outdir, job.name[1:])
      
      bleeder = Bleeder(jobdict, self._pools(len(kpoints), comm), comm)
      for value, job in bleeder.itercompute(): continue
      bleeder.cleanup()

      return Extract(outdir, comm)

  # otherwise, will drop out on call.
  else: 
    def __call__(self, *args, **kwargs):
      raise ImportError('Cannot use KEscan without jobs package.')
  
  def _interpret_kpoints(self, kpoints, vffout):
     """ Returns list of kpoints. """
     from numpy import zero, array
     # case where kpoints is None.
     if kpoints == None: return [zeros((3,1), dtype='float64')]
     # case where kpoints is already a single vector.
     if hasattr(kpoints, '__len__'):
       if len(kpoints) == 0: return [zeros((3,1), dtype='float64')]
       if len(kpoints) == 3:
         if not hasattr(kpoints[0], '__len__'): return [array(kpoints, dtype='float64')]
         if len(kpoints[0]) == 1: return [array([k[0] for k in kpoints], dtype='float64')]
     # case where kpoints is a callable.
     if hasattr(kpoints, '__call__'): kpoints = kpoints(vffout.structure)
     # last case covers list of vectors and finishes up callable.
     result = []  
     for k in kpoints: 
       kpoint = array(k, dtype='float64')
       assert len(kpoint) == 3, ValueError('k-vector = {0}?'.format(kpoint))
       result.append(kpoint)
     return result

  def _pools(self, N, comm):
    """ Optimizes number of pools. 
    
        Tries to find the largest number of pools which divides the number of
        kpoints, of procs, and the fft mesh, in that order increasingly
        inclusive order. Returns 1 on failure.
    """
    if comm == None: return 1
    if comm.size == 1: return 1
    if N == 1: return 1

    # finds divisors of N
    pools = [i for i in range(1, N+1) if N % i == 0 and comm.size / i > 0]
    if len(pools) == 1: return pools[-1]
    if len(pools) == 0: return 1

    # checks for pools which best divide comm.size
    pools = [i for i in pools if comm.size % i == 0]
    if len(pools) == 1: return pools[-1]
    if len(pools) == 0: return 1

    # checks for pools which best divide mesh size
    fftsize = self.fft_mesh[0] * self.fft_mesh[1] * self.fft_mesh[2]
    pools = [i for i in pools if fftmesh % (comm.size / i) == 0]
    if len(pools) >= 1: return pools[-1]
    return 1


class KGrid(KPoints):
  """ Unreduces kpoint grid with offsets. """

  def __init__(self, grid = None, offset = None):
    """ Initializes unreduced KGrid. """
    from numpy import array
    KPoints.__init__(self)
    self.grid = grid if grid != None else array([1,1,1])
    """ Grid dimensions in reciprocal space. """
    self.offset = offset if offset != None else array([0,0,0])
    """ Offset from Gamma of the grid. """

  def _mnk(self, structure):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv, norm
    from numpy import zeros, array, dot
    from ..crystal import into_voronoi
    inv_cell = structure.cell.T
    cell = inv(inv_cell)
    a = zeros((3,), dtype='float64')
    seen = []
    weight = 1e0 / float(self.grid[0] * self.grid[1] * self.grid[2])
    for x in xrange(self.grid[0]):
      a[0] = float(x + self.offset[0]) / float(self.grid[0])
      for y in xrange(self.grid[1]):
        a[1] = float(y + self.offset[1]) / float(self.grid[1]) 
        for z in xrange(self.grid[2]):
          a[2] = float(z + self.offset[2]) / float(self.grid[2])
          b = into_voronoi(dot(cell, a), cell, inv_cell)
          found = False
          for i, (count, v) in enumerate(seen):
            if all(abs(v-b) < 1e-12):
              seen[i][0] += weight
              found = True
              break
          if found == False: seen.append([weight, b.copy()])
    return seen
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_one = all( abs(array(self.grid)-array([1,1,1])) < 1e-12 )
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_one and is_zero:
      return '{0.__class__.__name__}()'.format(self)
    if is_one:
      return '{0.__class__.__name__}(offset=({0.offset[0]},{0.offset[0]},{0.offset[0]}))'\
             .format(self)
    if is_zero:
      return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}))'\
             .format(self)
    return '{0.__class__.__name__}(({0.grid[0]},{0.grid[0]},{0.grid[0]}), '\
           '({0.offset[0]},{0.offset[0]},{0.offset[0]}))'.format(self)


class KDensity(KPoints):
  """ Unreduced kpoint grid parameterized by the density and offset. """

  def __init__(self, density, offset = None):
    """ Initializes unreduced KGrid. """
    from numpy import array
    KPoints.__init__(self)
    self.density = density
    """ 1-dimensional density in cartesian coordinates (1/Angstrom). """
    self.offset = offset if offset != None else array([0,0,0])
    """ Offset from Gamma of the grid. """

  def _mnk(self, structure):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv, norm
    from numpy import zeros, array, dot, floor
    from ..crystal import fill_structure
    from ..crystal.gruber import Reduction
    from quantities import angstrom
    
    reduction = Reduction()
    cell = reduction(structure.cell, recip=True) * structure.scale
    density = self.density
    if hasattr(density, 'rescale'): density.rescale(1e0/Angstrom)
    grid = [0,0,0]

    for i in range(3):
      grid[i] = norm(cell[:,i]) / self.density
      grid[i] = int(max(1, floor(grid[i]+0.5)))

    kgrid = KGrid(grid=grid, offset=self.offset)
    return kgrid._mnk(structure)
 
  def __repr__(self):
    """ Represents this object. """
    from numpy import array, abs, all
    is_zero = all( abs(array(self.offset)-array([0,0,0])) < 1e-12 )
    if is_zero: return '{0.__class__.__name__}({0.density})'.format(self, repr(self.cell))
    return '{0.__class__.__name__}(({1}, ({2[0]},{2[0]},{2[0]}))'\
           .format(self, repr(self.cell), self.offset)

def _reduced_grids_factory(name, base):
  class ReducedKGrid(base): 
    """ Reduced {0} according to symmetries. """.format(base.__name__)
    def __init__(self, grid=None, offset=None, tolerance=1e-12):
      """ Initializes reduces k-grid. """
      base.__init__(self, grid, offset)
      self.tolerance = tolerance
      """ Criteria to determine whether two k-vectors are the same. """

    def _mnk(self, structure):
      """ Returns list of inequivalent vectors with multiplicity. """
      from numpy.linalg import inv, norm
      from ..crystal import Lattice, zero_centered
      lattice = Lattice()
      lattice.cell = inv(structure.cell.T)
      lattice.add_site = (0,0,0), '0'
      lattice.find_space_group()
      inv_cell = structure.cell.T
      # now checks whether symmetry kpoint exists or not.
      seen = []
      for mult, kpoint in base._mnk(self, structure): 
        found = False
        for i, (count, vec) in enumerate(seen):
          for op in lattice.space_group:
            u = zero_centered(vec-op(kpoint), lattice.cell, inv_cell)
            if all(abs(u)) < self.tolerance:
              found = True
              seen[i][0] += mult
              break
          if found: break
        if found == False: seen.append([mult, kpoint.copy()])
      return seen

    def __repr__(self):
      """ Represents this object. """
      if self.tolerance == 1e-12: return base.__repr__(self)
      result = base.__repr__(self)
      result = result[:-1].rstrip()
      if result[-1] == '(': return result + 'tolerance={0})'.format(self.tolerance)
      return result + ', tolerance={0})'.format(self.tolerance)
  ReducedKGrid.__name__ = name
  return ReducedKGrid


ReducedKGrid    = _reduced_grids_factory('ReducedKGrid', KGrid)
ReducedKDensity = _reduced_grids_factory('ReducedKDensity', KDensity)
