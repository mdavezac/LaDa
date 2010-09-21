""" Subpackage containing extraction methods for vasp parameters from vasp output. """
from ._dft import _ExtractImpl
from ._gw import _ExtractGWImpl

class Extract(_dft._ExtractImpl):
  """ Main class for extracting VASP output as python objects.

      This class should contain attributes (eg fermi_energy) which can extract
      their values from the vasp output files located in self.directory.  

      >>> result = Extract(directory = "./", comm = boost.mpi.world)
      >>> print result.fermi_energy * 13.26

      It would be preferable to limit these output files to L{files.OUTCAR}, as
      much as possible. Results are cached. To delete cache (and re-read results
      from output files), call C{self.uncache()}.
  """

  success         = property(_ExtractImpl._get_success)
  energy_sigma0   = property(_ExtractImpl._get_energy_sigma0)
  energy          = property(_ExtractImpl._get_energy)
  total_energy    = property(_ExtractImpl._get_energy)
  free_energy     = property(_ExtractImpl._get_free_energy)
  fermi_energy    = property(_ExtractImpl._get_fermi_energy)
  system          = property(_ExtractImpl._get_structure)
  structure       = property(_ExtractImpl._get_structure)
  contcar         = property(_ExtractImpl._get_contcar_structure)
  species         = property(_ExtractImpl._get_species)
  fft             = property(_ExtractImpl._get_fft)
  kpoints         = property(_ExtractImpl._get_kpoints)
  multiplicity    = property(_ExtractImpl._get_multiplicity)
  eigenvalues     = property(_ExtractImpl._get_eigenvalues)
  occupations     = property(_ExtractImpl._get_occupations)
  pressure        = property(_ExtractImpl._get_pressure)
  pulay_pressure  = property(_ExtractImpl._get_pulay_pressure)
  partial_charges = property(_ExtractImpl._get_partial_charges)
  magnetization   = property(_ExtractImpl._get_magnetization)
  moment          = property(_ExtractImpl._get_moment)
  nb_electrons    = property(_ExtractImpl._get_nb_electrons)
  ions_per_specie = property(_ExtractImpl._get_ions_per_type)
  name            = property(_ExtractImpl._get_name)
  vasp            = property(_ExtractImpl._get_functional)
  electropot      = property(_ExtractImpl._get_electropot)
  functional      = property(_ExtractImpl._get_functional)

  def __init__(self, directory = "", comm = None): 
    """ Initializes the extraction class. 

        @param comm: MPI group communicator. Extraction will be performed
                        for all procs in the group. In serial mode, comm can
                        be None.
        @param comm: boost.mpi.Communicator
        @param directory: path to the directory where the VASP output is located.
        @type directory: str
    """
    super(Extract, self).__init__(directory, comm)

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
    comm = self.comm 
    self.comm = None
    copy = deepcopy(self)
    self.comm = comm
    return copy

  def uncache(self): 
    """ Removes cached results.

        After this outputs are re-read from file.
    """
    from ...opt.decorators import uncache
    uncache(self)

class ExtractGW(_ExtractGWImpl):
  """ Extract data from VASP-GW output. """

  success         = property(_ExtractImpl._get_success)
  system          = property(_ExtractImpl._get_structure)
  structure       = property(_ExtractImpl._get_structure)
  species         = property(_ExtractImpl._get_species)
  kpoints         = property(_ExtractImpl._get_kpoints)
  multiplicity    = property(_ExtractImpl._get_multiplicity)
  dft_eigenvalues = property(_ExtractGWImpl._get_dft_eigenvalues)
  qp_eigenvalues  = property(_ExtractGWImpl._get_qp_eigenvalues)
  eigenvalues     = property(_ExtractGWImpl._get_qp_eigenvalues)
  self_energies   = property(_ExtractGWImpl._get_self_energies)
  occupations     = property(_ExtractGWImpl._get_occupations)

  def __init__(self, directory = "", comm = None): 
    """ Initializes the extraction class. 

        @param comm: MPI group communicator. Extraction will be performed
                        for all procs in the group. In serial mode, comm can
                        be None.
        @param comm: boost.mpi.Communicator
        @param directory: path to the directory where the VASP output is located.
        @type directory: str
    """
    super(ExtractGW, self).__init__(directory, comm)

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
    copy = deepcopy(self)
    return copy

try: from ... import jobs
except ImportError: pass
else: 
  class MassExtract(jobs.AbstractMassExtract):
    """ Propagates vasp extractors from all subdirectories.
    
        Trolls through all subdirectories for vasp calculations, and organises
        results as a dictionary where keys are the name of the diretory.
    """
    def __init__(self, path, hook = None, comm = None, Extract = None):
      """ Initializes MassExtract.
      
      
          :Parameters:
            - `path` : root directory for which to investigate all subdirectories.
            - `hook` : Optional callable taking a path as an argument and
              returning True or False. If it returns False, that particular
              directory and its subdirectories will not be listed as a vasp
              directory.
            - `comm` : boost.mpi.communicator. Can be None.
            - Extract : Extraction class to use. Defaults to `lada.vasp.Extract.`
      """
      from os.path import exists, isdir
      from . import Extract as VaspExtract
      self.hook = hook
      """ Optional callable to exclude directories from extraction. 
      
          Callable takes a path as an argument and returns True or False. If it
          returns False, that particular ill not be listed as a vasp directory.
      """
      self.Extract = Extract if Extract != None else VaspExtract
      """ Extraction class to use. """

      super(MassExtract, self).__init__(path, comm = None)
      assert exists(self.root), RuntimeError("Path {0} does not exist.".format(self.root))
      assert isdir(self.root), RuntimeError("Path {0} is not a directory.".format(self.root))

    def walk_through(self):
      """ Goes through all directories with a contcar. """
      from os import walk, getcwd
      from os.path import abspath, relpath, abspath, join

      for dirpath, dirnames, filenames in walk(self.root, topdown=True, followlinks=True):
        if "OUTCAR" not in filenames: continue
        if hasattr(self.hook, "__call__") and not self.hook(join(self.root, dirpath)):
          while len(dirnames): dirnames.pop()
          continue

        try: result = self.Extract(join(self.root, dirpath), comm = self.comm)
        except: pass
        else: yield relpath(dirpath, self.root), result

    def _properties(self): 
      """ Returns cached __dir__ result. """
      return set([u for u in dir(self.Extract()) if u[0] != '_'])

      


