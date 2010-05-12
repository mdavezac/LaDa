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
