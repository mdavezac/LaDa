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
  species         = property(_ExtractImpl._get_species)
  fft             = property(_ExtractImpl._get_fft)
  kpoints         = property(_ExtractImpl._get_kpoints)
  multiplicity    = property(_ExtractImpl._get_multiplicity)
  eigenvalues     = property(_ExtractImpl._get_eigenvalues)
  occupations     = property(_ExtractImpl._get_occupations)
  pressure        = property(_ExtractImpl._get_pressure)
  pulay_pressure  = property(_ExtractImpl._get_pulay_pressure)
  partial_charges = property(_ExtractImpl._get_partial_charges)
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
    comm = self.comm 
    self.comm = None
    copy = deepcopy(self)
    self.comm = comm
    return copy

