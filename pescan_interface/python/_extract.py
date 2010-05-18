""" Module to extract esca and vff ouput. """

from ..vff import Extract as VffExtract, _get_script_text
from ..opt.decorators import broadcast_result, make_cached

class Extract(object):
  """ A class to extract data from ESCAN output files. """
  def __init__(self, directory = None, comm = None, escan = None):
    """ Initializes ESCAN extraction class. """
    from . import Escan

    if escan == None: escan = Escan()
    super(Extract, self).__init__()
    
    self._vffout = VffExtract(directory, comm = comm, vff = escan.vff)
    """ Private reference to vff extraction object. """
    self.OUTCAR = escan.OUTCAR
    """ OUTCAR file to extract stuff from. """
    self.comm = comm
    """ Communicator for extracting stuff. 

        All procs will get same results at end of extraction. 
        Program will hang if not all procs are called when extracting some
        value. Instead, use L{solo}.

        >>> extract.success # Ok
        >>> if comm.rank == 0: extract.success # will hang if comm.size != 1
        >>> if comm.rank == 0: extract.solo().success # Ok
    """
    self.directory = directory if directory != None else getcwd()
  
  def _get_directory(self):
    """ Directory with VASP output files """
    return self._directory
  def _set_directory(self, dir):
    from os.path import abspath, expanduser
    dir = abspath(expanduser(dir))
    if hasattr(self, "_directory"): 
      if dir != self._directory: self.uncache()
    self._directory = dir
  directory = property(_get_directory, _set_directory)

  @property
  def structure(self): 
    """ Greps structure from self.L{vff}.L{OUTCAR} """
    return self._vffout.structure
  @property
  def vff(self): 
    """ Greps vff functional from self.L{vff}.L{OUTCAR} """
    return self._vffout.vff
  @property
  def lattice(self): 
    """ Greps lattice from self.L{vff}.L{OUTCAR} """
    return self._vffout.lattice
  @property
  def minimizer(self): 
    """ Greps minimizer from self.L{vff}.L{OUTCAR} """
    return self._vffout.minimizer
  @property
  def energy(self): 
    """ Greps strain energy from self.L{vff}.L{OUTCAR} """
    return self._vffout.energy
  @property
  def stress(self): 
    """ Greps stress from self.L{vff}.L{OUTCAR} """
    return self._vffout.stess

  @property
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks for VFF success.
        
        At this point, checks for files and 
    """
    from os.path import exists, join
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    if not exists(path): return False

    with open(path, "r") as file:
      good = 0
      for line in file:
        if line.find("FINAL eigen energies, in eV") != -1: good += 1
        if line.find("# Computed ESCAN in:") != -1 and good == 1: good += 1; break
    return good == 2

  @property
  @make_cached
  def escan(self):
    """ Greps escan functional from self.L{OUTCAR}. """
    from os.path import exists, join
    from . import Escan, localH, nonlocalH, soH, AtomicPotential
    from numpy import array
    
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (file))

    @broadcast_result(attr=True, which=0)
    def get_functional(this):
      with open(path, "r") as file: return _get_script_text(file, "Escan")
    local_dict = { "lattice": self.lattice, "minimizer": self.minimizer,\
                   "vff": self.vff, "Escan": Escan, "localH": localH,\
                   "nonlocalH": nonlocalH, "soH": soH, "AtomicPotential":AtomicPotential,\
                   "array": array }
    exec get_functional(self) in globals(), local_dict
    return local_dict["escan"]


  def _double_trouble(self):
    """ Returns true, if non-spin polarized or Kammer calculations. """
    from numpy.linalg import norm
    from . import soH
    potential = self.solo().escan.potential
    if potential != soH: return True
    return norm(self.solo().escan.kpoint) < 1e-12


  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def eigenvalues(self):
    """ Greps eigenvalues from self.L{OUTCAR}. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (file))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL eigen energies, in eV") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find("*********************************") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    return array(result, dtype="float64") if not self._double_trouble()\
           else array([result[i/2] for i in range(2*len(result))], dtype="float64")

  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def convergence(self):
    """ Greps eigenvalue convergence errors from self.L{OUTCAR}. 
    
        Always returns "spin-polarized" number of eigenvalues.
    """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (file))
    with open(path, "r") as file:
      for line in file: 
        if line.find(" FINAL err of each states, A.U") != -1: break
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")
      result = []
      for line in file:
        if line.find(" FINAL eigen energies, in eV") != -1: break
        result.extend( float(u) for u in line.split() )
      else: raise IOError("Unexpected end of file when grepping for eigenvectors.")

    return array(result, dtype="float64") if not self._double_trouble()\
           else array([result[i/2] for i in range(2*len(result))], dtype="float64")

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nnodes(self):
    """ Greps eigenvalue convergence errors from self.L{OUTCAR}. """
    from os.path import exists, join
    from numpy import array
    path = self.OUTCAR
    if len(self.directory): path = join(self.directory, self.OUTCAR)
    assert exists(path), RuntimeError("Could not find file %s:" % (file))
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
    if self._raw_gwfns_data == None: return None
    result = []
    if self.is_spinor:
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(i, eig, self.raw_gwfns[:,i/2,0],\
                                        self.raw_gwfns[:,i/2,1], attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(i, eig, self.raw_gwfns[self.inverse_indices,i/2,0],\
                                        self.raw_gwfns[self.inverse_indices,i/2,1], \
                                        attenuation = self.attenuation) )
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(i, eig, self.raw_gwfns[:,i,0],\
                                      self.raw_gwfns[:,i,1], attenuation = self.attenuation) )
    else: # no spin polarization.
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(i, eig, self.raw_gwfns[:,i/2,0],\
                                        attenuation = self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(i, eig, self.raw_gwfns[self.inverse_indices,i/2,0], \
                                        attenuation = self.attenuation) )
          result.append(result[-1])
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(i, eig, self.raw_gwfns[:,i,0],None, self.attenuation) )
          result.append(result[-1])
    return result

  @property
  @make_cached
  def rwfns(self):
    """ Creates list of rWavefuntion objects. """
    from ._wfns import rWavefunction, gtor_fourrier
    if self._raw_gwfns_data == None: return
    result = []
    if self.is_spinor:
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm), \
            gtor_fourrier( self.raw_gwfns[self.inverse_indices,:,:],\
                           self.rvectors, self.gvectors, self.comm )
        for i, eig in enumerate(self.eigenvalues):
          rwfn = rWavefunction(i, eig, self._raw_rwfns[i%2][:,i/2,0], self._raw_rwfns[i%2][:,i/2,1])
          result.append(rwfn)
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          rwfn = rWavefunction(i, eig, self._raw_rwfns[i%2][:,i,0], self._raw_rwfns[i%2][:,i,1])
          result.append(rwfn)
    else: # no spin polarization.
      if self.is_krammer:
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm), \
            gtor_fourrier( self.raw_gwfns[self.inverse_indices,:,:],\
                           self.rvectors, self.gvectors, self.comm )
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(i, eig, self._raw_rwfns[i%2][:,i/2,0]) )
      else: # no krammer degeneracy
        self._raw_rwfns = \
            gtor_fourrier(self.raw_gwfns, self.rvectors, self.gvectors, self.comm)
        for i, eig in enumerate(self.eigenvalues):
          result.append( rWavefunction(i, eig, self._raw_rwfns[:,i,0]) )
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
    
        This is a tuple described making up the return of
        L{read_wavefunctions<lada.escan._escan.read_wavefunctions>}
    """
    from os import remove
    from os.path import exists, join
    from numpy.linalg import norm
    from boost.mpi import world
    from ..opt import redirect
    from ..opt.changedir import Changedir
    from ._escan import read_wavefunctions
    from . import soH

    assert self.comm.size >= self.nnodes,\
           RuntimeError("Must read wavefunctions with at least "\
                        "as many nodes as they were written to disk.")
    if self.comm.size > self.nnodes:
      color = 0 if self.comm.rank < self.nnodes else 1
      local_comm = self.comm.split(color)
    else: color, local_comm = 0, self.comm
    if color == 1: return None
    with redirect(fout="") as streams:
      with Changedir(self.directory) as directory:
        assert exists(self.escan.WAVECAR),\
               IOError("%s does not exist." % (join(self.directory, self.escan.WAVECAR)))
        self.escan._write_incar(self.comm, self.structure)
        nbstates = self.escan.nbstates if self.escan.potential == soH and norm(self.escan.kpoint)\
                   else self.escan.nbstates / 2
        result = read_wavefunctions(self.escan, range(nbstates), local_comm)
        remove(self.escan._INCAR + "." + str(world.rank))
    return result

  @property
  def raw_gwfns(self):
    """ Raw wavefunction data in g-space. 
    
        Numpy array with three axis: (i) g-vectors, (ii) bands, (iii) spins:

        >>> self.raw_gwfns[:,0, 0] # spin up components of band index 0.
        >>> self.raw_gwfns[:,0, 1] # spin down components of band index 0.
        >>> for i, g in enumerate(self.gvectors): # looks for G=0 component
        >>>   if np.linalg.norm(g) < 1e-8: break
        >>> self.raw_gwfns[i, 0, 0] # G=0 component of spin-up wavefunction with band-index 0.

        The band index is the one from pescan, eg. it is different if user
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
  def _wavefunction_path(self): return self.solo().escan.WAVECAR

  @property
  def is_spinor(self):
    """ True if wavefunction is a spinor. """
    from . import soH
    return self.escan.potential == soH

  @property
  def is_krammer(self):
    """ True if wavefunction is a spinor. """
    from numpy.linalg import norm
    from . import soH
    return norm(self.escan.kpoint) < 1e-12


  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. C{self.L{solo}()} returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    from copy import deepcopy
    
    if self.comm == None: return self
    copy = Extract(self.directory, comm = None)
    copy.OUTCAR = self.OUTCAR
    copy._vffout = self._vffout.solo()
    return copy
