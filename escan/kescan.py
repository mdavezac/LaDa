""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ["KEscan", 'KPoints', 'KGrid', 'ReducedKGrid']
from abc import ABCMeta, abstractmethod
from .functional import Functional
from .. import __all__ as all_lada_packages
 
class KPoints(object):
  """ Abstract base class for callable KMesh objects. """
  __metaclass__ = ABCMeta
  @abstractmethod
  def __call__(self, structure):
    """ Generator yielding kpoints or function returning list of kpoints. """
    pass
  @abstractmethod
  def __repr__(self, structure):
    """ This object must be representable. """
    pass

class KEscan(Functional):
  """ A wrapper around Escan for computing many kpoints. """
  def __init__(self, kpoints = None, multiplicity=None):
    """ Initializes the KEscan functional. """
    self.kpoints = kpoints
    """ Kpoints to use for calculations.
    
        This object must be None (Gamma), a KMesh-derived structure, a single
        kpoint, or a list of kpoints.
    """

  # need jobs package to run this code.
  if 'jobs' in all_lada_packages: 
    def __call__(self.structure, outdir=None, comm=None, *args, **kwargs):
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
        if hasattr(j
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
    def __call__(*args, **kwargs):
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
    from numpy import array, ide
    self.grid = grid if grid != None else array([1,1,1])
    """ Grid dimensions in reciprocal space. """
    self.offset = offset if offset != None else array([0,0,0])

  def __call__(self, structure):
    """ Yields kpoints on the grid. """
    from numpy.linalg import inv
    from numpy import zeros
    cell = inv(structure.cell.T)
    a = zeros((3,), dtype='float64')
    for x in xrange(self.grid[0]):
      a[0] = float(x) / float(self.grid[0]) + self.offset[0]
      for y in xrange(self.grid[1]):
        a[1] = float(y) / float(self.grid[1]) + self.offset[1]
        for z in xrange(self.grid[2]):
          a[2] = float(z) / float(self.grid[2]) + self.offset[2]
          yield cell * a

class ReducedKGrid(KGrid): 
  """ Reduces KGrid according to symmetries. """
  def __init__(self, grid=None, offset=None, tolerance=1e-12):
    """ Initializes reduces k-grid. """
    KGrid.__init__(self, grid, offset)
    self.tolerance = tolerance
    """ Criteria to determine whether two k-vectors are the same. """

  def __call__(self, structure):
    """ Yields KPoint grid reduced by symmetry. """
    from numpy.linalg import inv
    from ..crystal import Lattice
    # creates reciprocal space lattice to get symmetries.
    lattice = Lattice()
    lattice.cell = inv(structure.cell.T)
    lattice.add_site = (0,0,0), '0'
    lattice.find_space_group()
    # now checks whether symmetry kpoint exists or not.
    seen = []
    for kpoint in self.kpoints: 
      if any(norm(s-op(kpoint)) < tolerance for s in seen for op in lattice.space_group):
        continue
      seen.append(kpoint)
      yield seen[-1]

  def multiplicity(self, structure):
    """ Yields KPoint grid reduced by symmetry. """
    from numpy.linalg import inv
    from ..crystal import Lattice
    lattice = Lattice()
    lattice.cell = inv(structure.cell.T)
    lattice.add_site = (0,0,0), '0'
    lattice.find_space_group()
    # now checks whether symmetry kpoint exists or not.
    seen = []
    for kpoint in self.kpoints: 
      found = False
      for i, (count, vec) in enumerate(seen):
        if any(norm(vec-op(kpoint)) < tolerance for op in lattice.space_group):
          seen[i] += 1
          found = True
          break
      seen.append((1, kpoint))
    return [count for count, vec in seen]
