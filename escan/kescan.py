""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ["KEscan"]
from abc import ABCMeta, abstractmethod
from .functional import Functional
from .. import __all__ as all_lada_packages
 
class KMesh(object):
  """ Abstract base class for callable KMesh objects. """
  __metaclass__ = ABCMeta
  @abstractmethod
  def __call__(self, vffrun = None):
    """ Generator yielding kpoints or function returning list of kpoints. """
    pass
  @abstractmethod
  def __repr__(self, structure):
    """ This object must be representable. """
    pass

class KEscan(Functional):
  """ A wrapper around Escan for computing many kpoints. """
  def __init__(self, kpoints = None):
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
      bleeder.cleanup()`

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
     if hasattr(kpoints, '__call__'): kpoints = kpoints(vffout)
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


