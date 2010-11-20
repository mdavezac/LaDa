""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ['KEscan', 'Extract']
from abc import ABCMeta, abstractmethod
from .functional import Escan
from .. import __all__ as all_lada_packages
from ..opt import AbstractExtractBase
from ._extract import MassExtract as EscanMassExtract
 
class KEscan(Escan):
  """ A wrapper around Escan for computing many kpoints. """
  def __init__(self, kpoints=None, multiplicity=None, **kwargs):
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
    from .kpoints import KContainer
    # case for simple containers.
    if kpoints == None: kpoints, multiplicity = [[0,0,0]], [1]
    if not hasattr(kpoints, '__call__'): self.kpoints = KContainer(kpoints, multiplicity)
    Escan.__init__(self)
    """ True if performing calculation on single point. """

  # need jobs package to run this code.
  if 'jobs' in all_lada_packages: 
    def __call__(self, structure, outdir=None, comm=None, _in_call=False, **kwargs):
      """ Performs calculcations. """
      if _in_call == True: # single calculation.
        return Escan.__call__(self, structure, outdir, comm, **kwargs)

      from os.path import join
      from ..jobs import JobDict, Bleeder

      kpoints = kwargs.pop('kpoints', self.kpoints)
      do_relax_kpoint = not getattr(kpoints, 'do_relax_kpoint', not self.do_relax_kpoint)
      do_relax_kpoint = kwargs.pop('do_relax_kpoint', do_relax_kpoint)
      is_mpi = False if comm == None else comm.size > 1
      is_root = True if not is_mpi else comm.rank == 0

      # performs vff calculations
      vffrun = kwargs.get('vffrun', None)
      if vffrun == None: 
        vffrun = Escan.__call__(self, structure, outdir, comm, do_escan=False, **kwargs)
        kwargs.pop('vffrun', None)
        if kwargs.get('genpotrun', None) == None: kwargs.pop('genpotrun', None)
  
      # create list of kpoints.
      kpoints = self._interpret_kpoints(kpoints, vffrun)
      # checks for 
      if len(kpoints) == 1:
        return self(structure, outdir, comm, _in_call=True, kpoint=kpoints[0], **kwargs)

      jobdict = JobDict()
      for i, kpoint in enumerate(kpoints):
        job = jobdict / 'kpoint_{0}'.format(i)
        job.functional = self
        job.jobparams['kpoint'] = kpoint
        job.jobparams['structure'] = structure
        job.jobparams['do_relax_kpoint'] = do_relax_kpoint
        job.jobparams['outdir'] = join(outdir, job.name[1:])
        job.jobparams['_in_call'] = True
        if kwargs.get('genpotrun', None) == None: job.jobparams['genpotrun'] = vffrun
        if kwargs.get('vffrun', None) == None:    job.jobparams['vffrun']    = vffrun
      
      bleeder = Bleeder(jobdict, self._pools(len(kpoints), comm), comm)
      for result, job in bleeder.itercompute(**kwargs): continue
      bleeder.cleanup()

      result = Extract(outdir, comm)
      result.jobdict = jobdict
      return result

  # otherwise, will drop out on call.
  else: 
    def __call__(self, *args, **kwargs):
      raise ImportError('Cannot use KEscan without jobs package.')
  
  def _interpret_kpoints(self, kpoints, vffout):
     """ Returns list of kpoints. """
     from numpy import zeros, array
     # case where kpoints is None.
     if kpoints == None: return [zeros((3,1), dtype='float64')]
     # case where kpoints is already a single vector.
     if hasattr(kpoints, '__len__'):
       if len(kpoints) == 0: return [zeros((3,1), dtype='float64')]
       if len(kpoints) == 3:
         if not hasattr(kpoints[0], '__len__'): return [array(kpoints, dtype='float64')]
         if len(kpoints[0]) == 1: return [array([k[0] for k in kpoints], dtype='float64')]
     # case where kpoints is a callable.
     if hasattr(kpoints, '__call__'):
       kpoints = kpoints.kpoint(vffout.input_structure, vffout.structure)
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

  def __repr__(self):
    """ Represents KEscan instance. """
    if not hasattr(self.kpoints, '__call__'): return Escan.__repr__(self)
    return 'from {0.kpoints.__class__.__module__} import {0.kpoints.__class__.__name__}\n'\
           .format(self) + Escan.__repr__(self)

class Extract(EscanMassExtract):
  """ Extraction class for KEscan. """

  def __init__(self, directory=None, comm=None, unreduce=True, **kwargs):
    """ Initializes the extraction object. """
    EscanMassExtract.__init__(self, directory, comm=comm, **kwargs)
    self.unreduce = unreduce
    """ Unreduced kpoints if True and if kpoints scheme sports a mapping method. """

  def _do_unreduce(self):
    """ True if should unreduce kpoints. """
    if self.unreduce == False: return False
    return hasattr(self.functional.kpoints, 'mapping')

  def __iter_alljobs__(self):
    """ Yields extraction objects for KEscan calculations. """
    from glob import iglob, glob
    from re import compile
    from os.path import isdir

    regex = compile(r'kpoint_(\d+)/')
    paths = [path for path in iglob('kpoint_*/') if isdir(path) and regex.search(path) != None]
    vffout = self.Extract()._vffout

    for path in sorted(paths, key=lambda x: int(regex.search(path).group(1))):
      inpaths = glob('*')
      dirnames =  [u for u in inpaths if isdir(u)]
      filenames = set(inpaths) - set(dirnames)
      if not self.__is_calc_dir__(path, dirnames, filename): continue

      try: result = self.Extract(join(self.rootdir, dirpath), comm = self.comm)
      except: continue

      result._vffout = vffout
      result.OUTCAR = self.OUTCAR
      result.FUNCCAR = self.FUNCCAR
      yield join('/', relpath(dirpath, self.rootdir)), result

  def __getitem__(self, name):
    """ Forks between integer and str keys. """
    if isinstance(name, int):
      name = 'escan_{0}'.format(self._index(name))
    return EscanMassself['escan_{0}'.format(value)]
 
  def __len__(self):
    """ Number of kpoint calculations. """
    if self._do_unreduce: return len([k for k in self.functional.mapping()])
    return len(self.items())

  def _index(self):
    """ Returns index, accounting for possible unreduce. """
    if name < 0: name += self.__len__()
    if name >= self.__len__() or name < 0:
      raise IndexError('Index out-of-range:{0}.'.format(name))
    return list(self.functional.mapping())[name] if self._do_unreduce else name

  def __iter__(self): 
    """ Iterates through individual kpoint calculations. """
    for i in range(len(self)): yield self[i]

  @property
  def kpoints(self):
    """ kpoint values. """
    from numpy import array
    return array((job.kpoint for job in self), dtype='float64')

  @property
  def multiplicity(self):
    """ Multiplicity of the kpoints. """
    from numpy import array, ones
    if self._do_unreduce: 
      return array((m for m in self.functional.kpoints.multiplicity), dtype='float64')
    else: return ones((len(self),), dtype='float64') / float(len(self))

  @property
  def eigenvalues(self):
    """ Eigenvalues across all kpoints. """
    from numpy import array
    if len(self) == 0: return array()
    return array((job.eigenvalues for job in self), dtype='float64') * self[0].eigenvalues.units

  @property
  def escan(self):
    """ Returns functional used for calculation. """
    for job in self.iteritems(): return job
  functional = escan

  @property
  def vbm(self): 
    """ Returns energy at vbm. """
    from numpy import array, max
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.vff.structure)
    units = self.eigenvalues.itervalues().next().units
    return max(array(self.eigenvalues.values())[:, nbe-2:nbe]) * units

  @property
  def cbm(self): 
    """ Returns energy at vbm. """
    from numpy import array, min
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.vff.structure)
    units = self.eigenvalues.itervalues().next().units
    return min(array(self.eigenvalues.values())[:, nbe:nbe+2]) * units

  @property 
  def directness(self):
    """ Difference in energy between the CBM at Gamma and the LUMO. """
    from numpy.linalg import norm
    from ..crystal import nb_valence_states
    lumo = self.cbm
    gamma = min((job for job in self.values()), key=lambda x: norm(x.escan.kpoint))
    if norm(gamma.escan.kpoint) > 1e-6: raise RuntimeError("Gamma point not found.")
    nbe = nb_valence_states(self.vff.structure)
    cbm = min(gamma.eigenvalues[nbe], gamma.eigenvalues[nbe+1])
    return cbm - lumo
