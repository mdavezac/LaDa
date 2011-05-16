""" Escan wrapper to compute many eigen k-points. """
__docformat__ = "restructuredtext en"
__all__ = ['KEscan', 'Extract']
from .functional import Escan
from .. import __all__ as all_lada_packages
from ..opt import AbstractExtractBase
from ..opt.decorators import make_cached
from ._extract import Extract as EscanExtract

class Extract(AbstractExtractBase):
  """ Extraction class for KEscan. """
  EscanExtract = staticmethod(EscanExtract)
  """ Escan extraction object. """

  def __init__(self, directory=None, comm=None, unreduce=True, **kwargs):
    """ Initializes the extraction object. """
    AbstractExtractBase.__init__(self, directory, comm=comm)
    self.unreduce = unreduce
    """ Unreduced kpoints if True and if kpoints scheme sports a mapping method. """

  @property
  def _do_unreduce(self):
    """ True if should unreduce kpoints. """
    if self.unreduce == False: return False
    return hasattr(self.functional.kpoints, 'mapping')

  @property
  def _joblist(self):
    """ List of cached jobs. 

        This property will return a reduced or unreduced list of cached jobs.
    """
    if '_cached_joblist' not in self.__dict__:
      from glob import iglob
      from re import compile
      from os.path import isdir, join, basename
      
      regex = compile(r'kpoint_(\d+)/')
      paths = [ path for path in iglob(join(self.directory, 'kpoint_*/'))\
                if isdir(path) and regex.search(path) != None ]
      paths = sorted(paths, key=lambda x: int(regex.search(x).group(1)))
      vffout = self.EscanExtract(self.directory, comm=self.comm)._vffout
      OUTCAR = self.EscanExtract().OUTCAR
      
      result = []
      for path in paths:
        filenames = [basename(u) for u in iglob(join(path, '*')) if not isdir(u)]
        if OUTCAR not in filenames: continue
      
        try: extractor = self.EscanExtract(path, comm = self.comm)
        except: continue
      
        extractor._vffout = vffout
        result.append(extractor)
      self.__dict__['_cached_joblist'] = result
    
    if self._do_unreduce:
      if '_cached_ujoblist' not in self.__dict__:
        result = []
        for i in self.functional.kpoints.mapping(self.input_structure, self.structure):
          assert i < len(self._cached_joblist),\
                 RuntimeError('{0}, {1}, {2}'.format(self, i, len(result)))
          result.append(self._cached_joblist[i])
        self.__dict__['_cached_ujoblist'] = result
      return self._cached_ujoblist

    return self._cached_joblist
        

  @property 
  def success(self):
    """ True if jobs are successfull. """
    try: 
      if len(self._joblist) != len(self.multiplicities): return False
      return all(job.success for job in self)
    except: return False

  @property
  def vff(self):
    """ Returns vff functional. """
    return self._rootrun.vff
  @property
  def structure(self):
    """ Returns vff output structure. """
    return self._rootrun.structure
  @property
  def input_structure(self):
    """ Returns input output structure. """
    return self._rootrun.input_structure
  @property
  def energy(self):
    """ Returns VFF energy. """
    return self._rootrun.energy
  @property
  def stress(self):
    """ Returns VFF energy. """
    return self._rootrun.stress

  def __getitem__(self, index):
    """ Returns extraction object of given kpoint calculation. """
    return self._joblist[index]
 
  def __len__(self):
    """ Number of kpoint calculations. """
    return len(self._joblist)

  def __iter__(self): 
    """ Iterates through individual kpoint calculations. """
    return self._joblist.__iter__()

  @property
  def _knm(self):
    """ kpoint and multiplicities. """
    if self._do_unreduce:
      if '_cached_uknm' not in self.__dict__:
        kpoints = self.functional.kpoints
        istr, ostr = self.input_structure, self.structure
        self.__dict__['_cached_uknm'] = [(m,k) for m, k in kpoints.unreduced(istr, ostr)]
      return self._cached_uknm
    else: 
      if '_cached_knm' not in self.__dict__:
        kpoints = self.functional.kpoints
        istr, ostr = self.input_structure, self.structure
        self.__dict__['_cached_knm'] = [(m,k) for m, k in kpoints(istr, ostr)]
      return self._cached_knm

  @property
  def kpoints(self):
    """ Kpoints in the escan run. 

        Depending on ``unreduce`` attribute, the kpoints are unreduced or not.
    """
    from numpy import array
    return array([k for m, k in self._knm])

  @property 
  def _computed_kpoints(self):
    """ Kpoints as grepped from calculation.

        These k-points include relaxation. They are given as computed in
        calculation, rather than from the functional. Probably more of a
        debugging value.
    """
    from numpy import array
    return array([u.kpoint for u in self])

  @property
  def multiplicity(self):
    """ Multiplicity of the kpoints. """
    from numpy import array
    return array([m for m, k in self._knm])
  multiplicities = multiplicity

  @property
  def eigenvalues(self):
    """ Eigenvalues across all kpoints. """
    from numpy import array
    from quantities import eV
    return array([job.eigenvalues.rescale(eV) for job in self], dtype='float64') * eV

  @property
  @make_cached
  def functional(self):
    """ Returns functional used for calculation. """
    return EscanExtract(self.directory, comm=self.comm).functional

  @property
  def vbm(self): 
    """ Energy at valence band minimum. """
    if self.functional.eref != None:
      raise RuntimeError('Cannot extract VBM from folded spectrum calculation.')
    from numpy import max
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.structure)
    return max(self.eigenvalues[:, nbe-2:nbe])

  @property
  def cbm(self): 
    """ Energy at conduction band minimum. """
    if self.functional.eref != None:
      raise RuntimeError('Cannot extract CBM from folded spectrum calculation.')
    from numpy import min
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.structure)
    return min(self.eigenvalues[:, nbe:nbe+2])

  @property
  def gap(self):
    """ Gap between the VBM and CBM. """
    return self.cbm - self.vbm

  @property
  def cbm_direct_gap(self):
    """ Gap between the CBM and valence band at the same kpoint. """
    if self.functional.eref != None:
      raise RuntimeError('Cannot extract CBM from folded spectrum calculation.')
    from numpy import argmin
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.structure)
    i = argmin(self.eigenvalues[:, nbe:nbe+2])
    return self.eigenvalues[:, nbe:nbe+2].flat[i] - self.eigenvalues[:,nbe-2:nbe].flat[i]

  @property
  def vbm_direct_gap(self):
    """ Gap between the VBM and conduction band at the same kpoint. """
    if self.functional.eref != None:
      raise RuntimeError('Cannot extract VBM from folded spectrum calculation.')
    from numpy import argmax
    from ..crystal import nb_valence_states
    nbe = nb_valence_states(self.structure)
    i = argmax(self.eigenvalues[:, nbe-2:nbe])
    return self.eigenvalues[:, nbe:nbe+2].flat[i] - self.eigenvalues[:,nbe-2:nbe].flat[i]

  @property
  def gamma_gap(self):
    """ Gap at Gamma if computed. """
    if self.functional.eref != None:
      raise RuntimeError('Cannot extract VBM from folded spectrum calculation.')
    from numpy import min, max
    from ..crystal import nb_valence_states

    for kindex, kpoint in enumerate(self.kpoints):
      if all(abs(kpoint) < 1e-12): break
    if not all(abs(kpoint) < 1e-12):
      raise RuntimeError('Could not find Gamma in computed kpoints.')

    nbe = nb_valence_states(self.structure)
    return min(self.eigenvalues[kindex, nbe:nbe+2]) - max(self.eigenvalues[kindex, nbe-2:nbe])

  @property 
  def directness(self):
    """ Difference in energy between the CBM at Gamma and the LUMO. """
    from numpy.linalg import norm
    from ..crystal import nb_valence_states
    lumo = self.cbm
    gamma = min((job for job in self.values()), key=lambda x: norm(x.escan.kpoint))
    if norm(gamma.escan.kpoint) > 1e-6: raise RuntimeError("Gamma point not found.")
    nbe = nb_valence_states(self.structure)
    cbm = min(gamma.eigenvalues[nbe], gamma.eigenvalues[nbe+1])
    return cbm - lumo

  @property 
  def _rootrun(self):
    """ Vff extraction object. """
    if '__rootrun' not in self.__dict__:
      self.__rootrun = self.EscanExtract(self.directory, comm = self.comm) 
    return self.__rootrun

  def iterfiles(self, **kwargs):
    """ Iterates over output/input files.

        Parameters are passed on to internal escan calculations.
    """
    for file in self._rootrun.iterfiles(**kwargs): yield file
    for kpoint_calc in self: 
      for file in kpoint_calc.iterfiles(**kwargs): yield file

  def uncache(self):
    """ Uncaches extracted values. """
    super(KExtract, self).uncache()
    for name in ['_cached_knm', '_cached_uknm', '_cached_joblist', '_cached_ujoblist']:
      self.__dict__.pop(name, None)

 
class KEscan(Escan):
  """ A wrapper around Escan for computing many kpoints. """
  Extract = Extract
  def __init__(self, kpoints=None, multiplicity=None, nbpools=-1, **kwargs):
    """ Initializes the KEscan functional. """
    self.kpoints = kpoints
    """ Kpoints to use for calculations.
    
        This object must be None (Gamma), a KMesh-derived structure, a single
        kpoint, or a list of kpoints.

        It is expected that the kpoints are expressed in cartesian coordinates
        of the reciprocal space. There is no 2|pi|. See `escan.KPoints` for more
        information.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    from .kpoints import KContainer
    # case for simple containers.
    if kpoints == None: kpoints, multiplicity = [[0,0,0]], [1]
    if hasattr(kpoints, '__call__'): self.kpoints = kpoints
    else: self.kpoints = KContainer(kpoints, multiplicity)
    escan_copy = kwargs.pop("escan", None) 

    self.nbpools = nbpools
    """ Number of processor pools (over kpoints). 

        If 0 or negative, then will try to determine it from number of kpoints,
        processors, and fft mesh.
    """
    super(KEscan, self).__init__(**kwargs)

    if escan_copy != None: # copy constructor from Escan instance. 
      from copy import deepcopy
      self.__dict__.update(deepcopy(escan_copy.__dict__))


  # need jobs package to run this code.
  if 'jobs' in all_lada_packages: 
    def __call__(self, structure, outdir=None, comm=None, **kwargs):
      """ Performs calculcations. """
      from inspect import getargspec
      from copy import deepcopy
      from os.path import join
      from ..jobs import JobDict, Bleeder, SuperCall
      from ..mpi import Communicator

      this = deepcopy(self)
      do_relax_kpoint = kwargs.pop('do_relax_kpoint', kwargs.pop('do_relax_kpoints', None))
      arglist = getargspec(Escan.__call__)[0]
      for key, value in kwargs.iteritems():
        if key in arglist: continue
        assert hasattr(this, key), TypeError("Unexpected keyword argument {0}.".format(key))
        setattr(this, key, value)
      if do_relax_kpoint != None: this.do_relax_kpoint = do_relax_kpoint
      comm = Communicator(comm, with_world=True)

      if not kwargs.get('overwrite', False): 
        extract = KEscan.Extract(directory=outdir, escan=this, comm=comm)
        if extract.success: return extract

      # performs vff calculations
      vffrun = kwargs.pop('vffrun', None)
      genpotrun = kwargs.pop('genpotrun', None)
      if vffrun == None or genpotrun == None: 
        kargs = kwargs.copy() # makes sure we don't include do_escan twice.
        kargs['do_escan'] = False
        out = super(KEscan, this).__call__( structure, outdir, comm, 
                                            vffrun = vffrun, genpotrun=genpotrun, **kargs )
        if vffrun    == None: vffrun    = out
        if genpotrun == None: genpotrun = out

  
      # create list of kpoints.
      kpoints = this._interpret_kpoints(this.kpoints, vffrun)

      jobdict = JobDict()
      for i, kpoint in enumerate(kpoints):
        job = jobdict / 'kpoint_{0}'.format(i)
        job.functional = SuperCall(KEscan, this)
        job.jobparams['kpoint']          = kpoint
        job.jobparams['structure']       = structure
        job.jobparams['do_relax_kpoint'] = False
        job.jobparams['outdir']          = join(outdir, job.name[1:])
        job.jobparams['vffrun']          = vffrun
        job.jobparams['genpotrun']       = genpotrun
      
      directory = '.' if outdir == None else outdir
      bleeder = Bleeder(jobdict, this._pools(len(kpoints), comm), comm, directory=directory)
      for result, job in bleeder.itercompute(**kwargs): continue
      bleeder.cleanup()

      result = Extract(outdir, comm, unreduce=True)
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
    from ..mpi import Communicator
    comm = Communicator(comm)
    if self.nbpools > 0: return min(comm.size, self.nbpools)
    if not comm.is_mpi: return 1
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
    pools = [i for i in pools if fftsize % (comm.size / i) == 0]
    if len(pools) >= 1: return pools[-1]
    return 1

  def __repr__(self):
    """ Represents KEscan instance. """
    if not hasattr(self.kpoints, '__call__'): return Escan.__repr__(self)
    return 'from {0.kpoints.__class__.__module__} import {0.kpoints.__class__.__name__}\n'\
           .format(self) + Escan.__repr__(self)

  @property
  def do_relax_kpoint(self):
    """ Whether to deform kpoints from original to relaxed geometry.
    
        Default is True. Relaxed cell is taken from `_POSCAR`
        Coding: Also sets attribute in kpoints. 
    """
    return self.__dict__["do_relax_kpoint"]
  @do_relax_kpoint.setter
  def do_relax_kpoint(self, value): 
    self.__dict__["do_relax_kpoint"] = value
    self.kpoints.relax = value

  do_relax_kpoints = do_relax_kpoint
  """ Alias for `KEscan.do_relax_kpoint`. """

  def to_escan(self):
    """ Converts to an `Escan` functional. """
    from .functional import Escan
    result = Escan()
    result.__dict__.update(self.__dict__)
    result.__dict__.pop('kpoints')
    return result
