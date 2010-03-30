""" Subpackage containing extraction methods for vasp parameters from vasp output. """
from ._dft import _ExtractImpl
from ._gw import _ExtractGWImpl
from ._success import Success

class Extract(_ExtractImpl):
  """ Main class for extracting VASP output as python objects.

      This class should contain attributes (eg fermi_energy) which can extract
      their values from the vasp output files located in self.directory.  

      >>> result = Extract(directory = "./", mpicomm = boost.mpi.world)
      >>> print result.fermi_energy * 13.26

      It would be preferable to limit these output files to L{files.OUTCAR}, as
      much as possible. Results are cached. To delete cache (and re-read results
      from output files), call C{self.uncache()}.
  """

  success = Success()
  r""" Checks for success of vasp calculation """
  energy_sigma0 = property(_ExtractImpl._get_energy_sigma0)
  r""" Gets total energy extrapolated to $\sigma=0$ from vasp run """
  energy = property(_ExtractImpl._get_energy)
  """ Gets total energy from vasp run 
     
      same as self.L{total_energy}
  """
  total_energy = property(_ExtractImpl._get_energy)
  """ Gets total energy from vasp run
     
      same as self.L{energy}
  """
  free_energy = property(_ExtractImpl._get_free_energy)
  r""" Gets total free energy from vasp run """
  fermi_energy = property(_ExtractImpl._get_fermi_energy)
  r""" Gets fermi energy from vasp run """
  system = property(_ExtractImpl._get_structure)
  r""" Returns the relaxed structure (L{lada.crystal.Structure}) as obtained from the CONTCAR. """
  structure = property(_ExtractImpl._get_structure)
  r""" Alias for self.L{system}. """
  species = property(_ExtractImpl._get_species)
  """ Atomic species in system. """
  fft = property(_ExtractImpl._get_fft)
  r""" Gets recommended fft grid (for wavefunctions). """
  kpoints = property(_ExtractImpl._get_kpoints)
  r""" A mx3 matrix where each row corresponds to a kpoint.
  
       kpoints are in cartesian units.
   """
  multiplicity = property(_ExtractImpl._get_multiplicity)
  r""" A mx1 matrix where each correspond to a k-point multiplicity """
  eigenvalues = property(_ExtractImpl._get_eigenvalues)
  r""" A matrix of eigenvalues where each row contains eigenvalues of a single kpoint. """
  occupations = property(_ExtractImpl._get_occupations)
  r""" A matrix of eigenvalues where each row contains band occupations of a single kpoint. """
  pressure = property(_ExtractImpl._get_pressure)
  """ External pressure at end of calculation. """
  pulay_pressure = property(_ExtractImpl._get_pulay_pressure)
  """ Pulay pressure at end of calculation. """
  partial_charges = property(_ExtractImpl._get_partial_charges)
  """ Partial charges.

      This is a numpy array where the first dimension is the ion (eg one row
      per ion), and the second the partial charges for each angular momentum.
      The total is not included.
  """

  def __init__(self, directory = "", mpicomm = None): 
    """ Initializes the extraction class. 

        @param mpicomm: MPI group communicator. Extraction will be performed
                        for all procs in the group. In serial mode, mpicomm can
                        be None.
        @param mpicomm: boost.mpi.Communicator
        @param directory: path to the directory where the VASP output is located.
        @type directory: str
    """
    super(Extract, self).__init__(directory, mpicomm)

  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. C{self.L{solo}()} returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    if self.mpicomm == None: return self
    return Extract(directory = self.directory, mpicomm = None)

  def uncache(self): 
    """ Removes cached results.

        After this outputs are re-read from file.
    """
    from _mpi import _uncache
    _uncache(self)

class ExtractGW(_ExtractGWImpl):
  """ Extract data from VASP-GW output. """

  success = Success()
  r""" Checks for success of vasp calculation """
  system = property(_ExtractImpl._get_structure)
  r""" Returns the relaxed structure (L{lada.crystal.Structure}) as obtained from the CONTCAR. """
  structure = property(_ExtractImpl._get_structure)
  r""" Alias for self.L{system}. """
  species = property(_ExtractImpl._get_species)
  """ Atomic species in system. """
  kpoints = property(_ExtractImpl._get_kpoints)
  r""" A mx3 matrix where each row corresponds to a kpoint.
  
       kpoints are in cartesian units.
   """
  multiplicity = property(_ExtractImpl._get_multiplicity)
  r""" A mx1 matrix where each correspond to a k-point multiplicity """
  dft_eigenvalues = property(_ExtractGWImpl._get_dft_eigenvalues)
  r""" A matrix of DFT eigenvalues.
  
       Each row corresponds to a k-point and each column to a band.
   """
  qp_eigenvalues = property(_ExtractGWImpl._get_qp_eigenvalues)
  r""" A matrix of qasi-particle eigenvalues.
  
       Each row corresponds to a k-point and each column to a band.
   """
  eigenvalues = property(_ExtractGWImpl._get_qp_eigenvalues)
  r""" Alias for L{qp_eigenvalues}. """
  self_energies = property(_ExtractGWImpl._get_self_energies)
  r""" A matrix of qasi-particle self-energies.
  
       Each row corresponds to a k-point and each column to a band.
   """
  occupations = property(_ExtractGWImpl._get_occupations)
  r""" A matrix of occupations.
  
       Each row corresponds to a k-point and each column to a band.
   """

  def __init__(self, directory = "", mpicomm = None): 
    """ Initializes the extraction class. 

        @param mpicomm: MPI group communicator. Extraction will be performed
                        for all procs in the group. In serial mode, mpicomm can
                        be None.
        @param mpicomm: boost.mpi.Communicator
        @param directory: path to the directory where the VASP output is located.
        @type directory: str
    """
    super(ExtractGW, self).__init__(directory, mpicomm)

  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. C{self.L{solo}()} returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    if self.mpicomm == None: return self
    return ExtractGW(directory = self.directory, mpicomm = None)

  def uncache(self): 
    """ Removes cached results.

        After this outputs are re-read from file.
    """
    from _mpi import _uncache
    _uncache(self)


