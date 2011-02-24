""" Module to extract esca and vff ouput. """
__docformat__ = "restructuredtext en"
__all__ = ['Extract', 'MassExtract']

from ..opt.decorators import broadcast_result, make_cached
from ..opt import AbstractExtractBase, OutcarSearchMixin
from ..jobs import AbstractMassExtractDirectories


class Extract(AbstractExtractBase, OutcarSearchMixin):
  """ A class to extract data from ESCAN output files. 
  
      This class helps to extract information from the escan output, including
      escan and vff parameters, relaxed crystal structure, eigenvalues, and
      optionally, wavefunctions in real and reciprocal space. Where possible
      quantities have attached units, using the package ``quantities``.
  """
  def __init__(self, directory = None, comm = None, escan = None):
    """ Initializes ESCAN extraction class. 
    
        :Parameters: 
          directory : str or None
            Directory where escan output files are located. Defaults to current
            working directory if None.
          comm : boost.mpi.communicator
            Communicator containing as many processes as were used to perform
            calculations. This is only mandatory when using wavefunctions in
            some way.
          escan : lada.escan.Escan
            Wrapper around the escan functional.
    """
    from ..vff import Extract as VffExtract
    from . import Escan

    super(Extract, self).__init__(directory=directory, comm=None)

    if escan == None: escan = Escan()
    
    self._vffout = VffExtract(directory, comm = None, vff = escan.vff)
    """ Private reference to vff extraction object. """
    self.OUTCAR = escan.OUTCAR
    """ OUTCAR file to extract stuff from. """
    self.FUNCCAR = escan._FUNCCAR
    """ Pickle to FUNCCAR. """
    self.comm = comm

  def __funccar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.FUNCCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  @property 
  def comm(self):
    """ Communicator over which to sync output. """
    return self._vffout.comm 
  @comm.setter
  def comm(self, value):
    if hasattr(self, "_vffout"): self._vffout.comm = value

  @property
  def kpoint(self):
    """ K-point in this calculation. """
    from numpy import array, zeros, pi
    from quantities import angstrom
    from ..physics import a0
    regex = r'\s*ikpt,akx,aky,akz\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)'
    result = self._find_first_OUTCAR(regex)
    assert result != None,\
           RuntimeError('Could not find kpoint in file {0};'.format(self.__outcar__().name))
    if result.group(1) == '0': return zeros((3,), dtype='float64')
    return array([result.group(2), result.group(3), result.group(4)], dtype='float64') / 2e0 / pi\
           * (self.structure.scale * angstrom).rescale(a0).magnitude

  
  def __directory__hook__(self):
    """ Called whenever the directory changes. """
    super(Extract, self).__directory_hook__()
    self._vffout.directory = value

  def uncache(self): 
    """ Uncache values. """
    super(Extract, self).uncache(self)
    self._vffout.uncache()

  def __copy__(self):
    """ Returns a shallow copy of this object. """
    result = super(Extract, self).__copy__()
    result._vffout = self._vffout.__copy__()
    return result

  def copy(self, **kwargs):
    """ Returns a shallow copy of this object.

        :param kwargs:
          Any keyword argument is set as an attribute of this object.
    """
    result = self.__copy__()
    for k, v in kwargs.items():
      if hasattr(result._vffout, k):
        setattr(result._vffout, k, v)
        if hasattr(result, k): setattr(result, k, v)
      else: setattr(result, k, v)
    return result

  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for Escan success.
        
        At this point, checks for files and 
    """
    from re import compile
    from os.path import exists, join
    do_escan_re = compile(r'functional\.do_escan\s*=')

    try:
      good = 0
      is_do_escan = True
      with self.__outcar__() as file:
        for line in file:
          if line.find("FINAL eigen energies, in eV") != -1: good += 1
          elif do_escan_re.search(line) != None: is_do_escan = eval(line.split()[-1])
          elif line.find("# Computed ESCAN in:") != -1: good += 1; break
      return (good == 2 and is_do_escan) or (good == 1 and not is_do_escan)
    except: return False

  @property
  @make_cached
  def functional(self):
    """ Greps escan functional from OUTCAR. """
    from os.path import exists, join
    from numpy import array
    from cPickle import load
    from ..opt.changedir import Changedir
    from ..vff import _get_script_text
    from . import Escan, localH, nonlocalH, soH, AtomicPotential
    
    # tries to read from pickle.
    path = self.FUNCCAR
    try:
      with self.__funccar__() as file: result = load(file)
    except: pass 
    else: return result


    # tries to read from outcar.
    @broadcast_result(attr=True, which=0)
    def get_functional(this):
      with self.__outcar__() as file: return _get_script_text(file, "Escan")
    local_dict = { "lattice": self.lattice, "minimizer": self.minimizer,\
                   "vff_functional": self.vff, "Escan": Escan, "localH": localH,\
                   "nonlocalH": nonlocalH, "soH": soH, "AtomicPotential":AtomicPotential,\
                   "array": array }
    # moves to output directory to get relative paths right.
    with Changedir(self.directory, comm=self.comm) as cwd:
      exec get_functional(self) in globals(), local_dict
    return local_dict["escan_functional"] if "escan_functional" in local_dict\
           else local_dict["functional"]

  @property 
  def escan(self): 
    """ Alias for functional. """
    from warnings import warn
    warn(DeprecationWarning('escan attribute is deprecated in favor of functional.'), stacklevel=2)
    return self.functional

  @property
  def _double_trouble(self):
    """ Returns true, if non-spin polarized or Kammer calculations. """
    from numpy.linalg import norm
    from . import soH
    seul = self.solo()
    if seul.functional.nbstates  ==   1: return False
    if seul.functional.potential != soH: return True
    return norm(seul.functional.kpoint) < 1e-12


  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def eigenvalues(self):
    """ Greps eigenvalues from OUTCAR. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    from quantities import eV
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file {0}.".format(path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL eigen energies, in eV") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find("*********************************") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    if self._double_trouble: result = [result[i/2] for i in range(2*len(result))]
    return array(result, dtype="float64") * eV

  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def convergence(self):
    """ Greps eigenvalue convergence errors from OUTCAR. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL err of each states, A.U") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find(" FINAL eigen energies, in eV") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    if self._double_trouble: result = [result[i/2] for i in range(2*len(result))]
    return array(result, dtype="float64") 

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nnodes(self):
    """ Greps eigenvalue convergence errors from OUTCAR. """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (path))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" nnodes =") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
    return int(line.split()[2])


  @property
  @make_cached
  def gwfns(self):
    """ Creates list of Wavefuntion objects. """
    from _wfns import Wavefunction
    result = []
    if self.is_spinor:
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i/2,0],\
                                        self.raw_gwfns[:,i/2,1], attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(self.comm, i, eig,\
                                        -self.raw_gwfns[self.inverse_indices,i/2,1].conjugate(),\
                                         self.raw_gwfns[self.inverse_indices,i/2,0].conjugate(), \
                                        attenuation = self.attenuation) )
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i,0],\
                                      self.raw_gwfns[:,i,1], attenuation = self.attenuation) )
    else: # no spin polarization.
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i/2,0],\
                                        attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(self.comm, i, eig, \
                                        self.raw_gwfns[self.inverse_indices,i/2,0], \
                                        attenuation = self.attenuation) )
          result.append(result[-1])
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(self.comm, i, eig, self.raw_gwfns[:,i,0],\
                                      None, self.attenuation) )
          result.append(result[-1])
    return result

  @property
  @make_cached
  def rwfns(self):
    """ Creates list of rWavefuntion objects. """
    from ._wfns import rWavefunction, gtor_fourrier
    result = []
    if self.is_spinor:
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          if i%2 == 0:
            rwfn = rWavefunction( self.comm, i, eig, self._raw_rwfns[:,i/2,0],\
                                  self._raw_rwfns[:,i/2,1])
          else: 
            rwfn = rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i/2,1].conjugate(),\
                                 -self._raw_rwfns[:,i/2,0].conjugate())
          result.append(rwfn)
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          rwfn = rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i,0], self._raw_rwfns[:,i,1])
          result.append(rwfn)
    else: # no spin polarization.
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          if i%2 == 0: 
            result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i/2,0]) )
          else:
            result.append( rWavefunction( self.comm, i, eig,
                                          -self._raw_rwfns[:,i/2,0].conjugate()) )
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(self.comm, i, eig, self._raw_rwfns[:,i,0]) )
    return result

  @property
  def raw_rwfns(self):
    """ Raw real-space wavefunction data. """
    if not hasattr(self, "_raw_rwfns"): self.rwfns # creates data
    return self._raw_rwfns

  @property
  @make_cached
  def _raw_gwfns_data(self):
    """ Reads and caches g-space wavefunction data. 
    
        This property is a tuple holding information about the wavefunctions.
        
        - a spin by N by x matrix holding the N wavefuntions/spinor.
        - a 3 by x matrix with each row a G-vector in units of
          `lada.physics.reduced_reciprocal_au`.
        - a 3 by x matrix with each row a R-vector in atomic units.
        - one-dimensional array of real coefficients to smooth higher energy G-vectors.
        - one-dimensional array of integer indices to map G-vectors to -G.
    """
    # check for mpi first
    from .. import lada_with_mpi
    assert lada_with_mpi, RuntimeError("Lada loaded without mpi. Cannot read wavefunctions.")
    # then check for function.
    from os.path import exists
    from numpy import sqrt
    from numpy.linalg import norm, det
    from quantities import angstrom, pi
    from boost.mpi import world
    from ..opt import redirect
    from ..opt.changedir import Changedir
    from ..physics import a0, reduced_reciprocal_au
    from ._escan import read_wavefunctions
    from . import soH

    comm = self.comm if self.comm != None else world
    is_root = comm.rank == 0 
    assert self.success
    assert self.nnodes == comm.size, \
           RuntimeError( "Must read wavefunctions with as many nodes as "\
                         "they were written to disk.")
    comm.barrier()
    with Changedir(self.directory, comm=self.comm) as directory:
      if is_root: 
        assert exists(self.functional.WAVECAR),\
               IOError("{0} does not exist.".format(self.functional.WAVECAR))
        if exists(self.functional._INCAR):
          from shutil import copyfile
          copyfile(self.functional._INCAR, ".lada_wavefunction_save_incar")
      self.functional._write_incar(self.comm, self.structure)
      if self.functional.potential == soH and norm(self.functional.kpoint):
        nbstates = self.functional.nbstates
      else: nbstates = self.functional.nbstates / 2
      with redirect(fout="") as streams:
        result = read_wavefunctions(self.functional, range(nbstates), comm, self.is_krammer)
      if is_root and exists(".lada_wavefunction_save_incar"):
        from shutil import move
        move(".lada_wavefunction_save_incar", self.functional._INCAR)
      elif is_root:
        from os import remove
        remove(self.functional._INCAR)
    comm.barrier()

    cell = self.structure.cell * self.structure.scale * angstrom
    normalization = det(cell.rescale(a0)) 
    return result[0] * sqrt(normalization), result[1] * 0.5 / pi * reduced_reciprocal_au,\
           result[2] * a0, result[3], result[4]

  @property
  def raw_gwfns(self):
    """ Raw wavefunction data in g-space. 
    
        Numpy array with three axis: (i) g-vectors, (ii) bands, (iii) spins:

        >>> self.raw_gwfns[:,0, 0] # spin up components of band index 0.
        >>> self.raw_gwfns[:,0, 1] # spin down components of band index 0.
        >>> for i, g in enumerate(self.gvectors): # looks for G=0 component
        >>>   if np.linalg.norm(g) < 1e-8: break
        >>> self.raw_gwfns[i, 0, 0] # G=0 component of spin-up wavefunction with band-index 0.

        The band index is the one from ESCAN, eg. it is different if user
        Krammer doubling or not, etc. This data is exactly as read from disk.
    """
    return self._raw_gwfns_data[0]

  @property
  def gvectors(self):
    """ G-vector values of wavefuntions. """
    return self._raw_gwfns_data[1]

  @property
  def rvectors(self):
    """ R-vector values of wavefuntions. """
    return self._raw_gwfns_data[2]

  @property
  def attenuation(self):
    """ G-vector attenuation values of wavefuntions. """
    return self._raw_gwfns_data[3]

  @property
  def inverse_indices(self):
    """ Indices to -G vectors of wavefuntions. """
    return self._raw_gwfns_data[4]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _wavefunction_path(self): return self.solo().functional.WAVECAR

  @property
  def is_spinor(self):
    """ True if wavefunction is a spinor. """
    from . import soH
    return self.functional.potential == soH

  @property
  def is_krammer(self):
    """ True if wavefunction has krammer degenerate equivalent. """
    return self.functional.is_krammer


  def __getattr__(self, name):
    """ Passes on public attributes to vff extractor, then to escan functional. """
    if name[0] != '_':
      if name in dir(self._vffout): return getattr(self._vffout, name)
      elif self.success and hasattr(self.functional, name):
        return getattr(self.functional, name)
    raise AttributeError("Unknown attribute {0}. It could be the "\
                         "run is unsuccessfull.".format(name))

  def __dir__(self):
    """ Returns list attributes.
    
        Since __getattr__ is modified, we need to make sure __dir__ returns a
        complete list of attributes. This is usefull  for command-line
        completion in ipython.
    """
    exclude = set(["add_potential", "write_escan_input"])
    result = [u for u in self.__dict__.keys() if u[0] != "_"]
    result.extend( [u for u in dir(self.__class__) if u[0] != "_"] )
    result.extend( [u for u in dir(self._vffout) if u[0] != "_"] )
    if self.success: result.extend( [u for u in dir(self.functional) if u[0] != "_"] )
    return list( set(result) - exclude )


class MassExtract(AbstractMassExtractDirectories):
  """ Extracts all escan calculations nested within a given input directory. """
  Extract = staticmethod(Extract)
  """ Extraction object for a single calculation. """
  def __init__(self, path = ".", **kwargs):
    """ Initializes AbstractMassExtractDirectories.
    
    
        :Parameters:
          path : str or None
            Root directory for which to investigate all subdirectories.
            If None, uses current working directory.
          kwargs : dict
            Keyword parameters passed on to AbstractMassExtractDirectories.

        :kwarg naked_end: True if should return value rather than dict when only one item.
        :kwarg unix_re: converts regex patterns from unix-like expression.
    """
    # this will throw on unknown kwargs arguments.
    if 'Extract' not in kwargs: kwargs['Extract'] = Extract
    super(MassExtract, self).__init__(path, **kwargs)
    del self.__dict__['Extract']

  def __iter_alljobs__(self):
    """ Goes through all directories with a contcar. """
    from os import walk, getcwd
    from os.path import abspath, relpath, abspath, join

    OUTCAR = self.Extract().OUTCAR
    for dirpath, dirnames, filenames in walk(self.rootdir, topdown=True, followlinks=True):
      if OUTCAR not in filenames: continue

      try: result = self.Extract(join(self.rootdir, dirpath), comm = self.comm)
      except: continue

      yield join('/', relpath(dirpath, self.rootdir)), result
      
  def __is_calc_dir__(self, dirpath, dirnames, filenames):
    """ Returns true this directory contains a calculation. """
    OUTCAR = self.Extract().OUTCAR
    return OUTCAR in filenames
