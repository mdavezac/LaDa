""" Subpackage containing extraction methods for vasp parameters from vasp output. """
from _success import Success
from _mpi import _bound_mpi_extraction
class Extract(object):
  """ Main class for extracting VASP output as python objects.

      This class should contain attributes (eg fermi_energy) which can extract
      their values from the vasp output files located in self.directory.  

      >>> result = Extract(directory = "./", mpicomm = boost.mpi.world)
      >>> print result.fermi_energy * 13.26

      It would be preferable to limit these output files to L{files.OUTCAR}, as
      much as possible. Results are cached. To delete cache (and re-read results
      from output files), call C{self.uncache()}.
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
    self.directory = directory
    """ path to the directory where the VASP output is located. """
    self.mpicomm = mpicomm
    """ MPI group communicator. """

  success = Success()
  r""" Checks for success of vasp calculation """

  def uncache(self): 
    from _mpi import _uncache
    _uncache(self)

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
        

  @_bound_mpi_extraction
  def _get_energy_sigma0(self):
    """ Gets total energy extrapolated to $\sigma=0$ from vasp run """
    from os.path import exists, join
    from re import compile, X as re_X
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                           energy\(sigma->0\)\s+=\s+(\S+)""", re_X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(2))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  energy_sigma0 = property(_get_energy_sigma0)
  r""" Gets total energy extrapolated to $\sigma=0$ from vasp run """

  @_bound_mpi_extraction
  def _get_energy(self):
    """ Gets total energy from vasp run """
    from os.path import exists, join
    from re import compile, X as re_X
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                           energy\(sigma->0\)\s+=\s+(\S+)""", re_X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  energy = property(_get_energy)
  """ Gets total energy from vasp run 
     
      same as self.L{total_energy}
  """
  total_energy = property(_get_energy)
  """ Gets total energy from vasp run
     
      same as self.L{energy}
  """

  @_bound_mpi_extraction
  def _get_free_energy(self):
    """ Gets total free energy from vasp run """
    from os.path import exists, join
    from re import compile
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""free\s+energy\s+TOTEN\s*=\s*(\S+)\s+eV""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  free_energy = property(_get_free_energy)
  r""" Gets total free energy from vasp run """

  @_bound_mpi_extraction
  def _get_fermi_energy(self):
    """ Gets total free energy from vasp run """
    from os.path import exists, join
    from re import compile
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""E-fermi\s*:\s*(\S+)""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result
  fermi_energy = property(_get_fermi_energy)
  r""" Gets fermi energy from vasp run """

  @_bound_mpi_extraction
  def _get_structure(self):
    """ Gets structure from CONTCAR file and total energy from OUTCAR """
    from os.path import exists, join
    from ..files import CONTCAR
    from ...crystal import read_poscar

    species_in = self.species

    path = CONTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)
    result = read_poscar(species_in, path)
    result.energy = self.energy
    return result
  system = property(_get_structure)
  r""" Returns the relaxed structure (L{lada.crystal.Structure}) as obtained from the CONTCAR. """
  structure = property(_get_structure)
  r""" Alias for self.L{system}. """

  @_bound_mpi_extraction
  def _get_species(self):
    """ Gets species from L{files.OUTCAR}. """
    from os.path import exists, join
    from re import compile, X as re_X
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      regex = compile(r"""VRHFIN\s*=\s*(\S+)\s*:""")
      for line in file:
        match = regex.search(line)
        if match == None: continue
        assert match.group(1) not in result, "Found the same element twice.\n" 
        result.append( match.group(1) )
    return tuple(result)
  species = property(_get_species)
  """ Atomic species in system. """

  @_bound_mpi_extraction
  def _get_fft(self):
    """ Returns recommended or actual fft setting """
    from os.path import exists, join
    from re import compile, search, X as re_X
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:

      # find start
      for line in file:
        if search("I would recommend the setting", line): break;
      ng_regex = compile(r"""WARNING:\s+wrap\s+around\s+error\s+
                                must\s+be\s+expected\s+set\s+NG(X|Y|Z)\s+to\s+(\d+)""", re_X)
      g_regex = compile(r"""NG(X|Y|Z)\s+is\s+ok\s+and\s+might\s+be\s+reduce\s+to\s+(\d+)""", re_X)
      found_regex = compile(r"""dimension\s+x,y,z\s+
                                   NGX\s+=\s+(\d+)\s+
                                   NGY\s+=\s+(\d+)\s+
                                   NGZ\s+=\s+(\d+)""", re_X)

      allset = 0
      fft = [None, None, None]
      for line in file:
        p = ng_regex.search(line)
        if p != None:
          if p.group(1) == 'X':
            fft[0] = int(p.group(2)) 
            allset += 1
          elif p.group(1) == 'Y':
            fft[1] = int(p.group(2))
            allset += 1
          elif p.group(1) == 'Z':
            fft[2] = int(p.group(2))
            allset += 1
          if allset == 3: break;
          continue;
        p = g_regex.search(line)
        if p != None:
          if p.group(1) == 'X':
            fft[0] = int(p.group(2)) 
            allset += 1
          elif p.group(1) == 'Y':
            fft[1] = int(p.group(2))
            allset += 1
          elif p.group(1) == 'Z':
            fft[2] = int(p.group(2))
            allset += 1
          if allset == 3: break;
          continue
        p = found_regex.search(line)
        if p != None:
          fft = [ int(p.group(1)), int(p.group(2)), int(p.group(3)) ]
          break;

      assert fft[0] != None, "File %s is incomplete or incoherent.\n" % (path)
      assert fft[1] != None, "File %s is incomplete or incoherent.\n" % (path)
      assert fft[2] != None, "File %s is incomplete or incoherent.\n" % (path)

      multiple = 8
      for i in range(3):
        if fft[i] % multiple: fft[i] += multiple - fft[i] % multiple
      return tuple(fft)
    raise RuntimeError, "File %s could not be opened.\n" % (path)
  fft = property(_get_fft)
  r""" Gets recommended fft grid (for wavefunctions). """


      
  @_bound_mpi_extraction
  def _get_kpoints(self):
    """ Returns kpoints. """
    from os.path import exists, join
    from re import compile, search 
    from numpy import array
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      found = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
      for line in file:
        if found.search(line) != None: break
      found = compile(r"""Following\s+cartesian\s+coordinates:""")
      for line in file:
        if found.search(line) != None: break
      file.next()
      for line in file:
        data = line.split()
        if len(data) != 4: break;
        result.append( data[:3] )
    return array(result, dtype="float64")
  kpoints = property(_get_kpoints)
  r""" A mx3 matrix where each row corresponds to a kpoint.
  
       kpoints are in cartesian units.
   """

  @_bound_mpi_extraction
  def _get_multiplicity(self):
    """ Returns multiplicity """
    from os.path import exists, join
    from re import compile, search 
    from ..files import OUTCAR
    from numpy import array

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      found = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
      for line in file:
        if found.search(line) != None: break
      found = compile(r"""Following\s+cartesian\s+coordinates:""")
      for line in file:
        if found.search(line) != None: break
      file.next()
      for line in file:
        data = line.split()
        if len(data) != 4: break;
        result.append( float(data[3]) )
    return array(result, dtype="float64")
  multiplicity = property(_get_multiplicity)
  r""" A mx1 matrix where each correspond to a k-point multiplicity """

  def _get_eigocc(self,which):
    """ Implementation of _get_eigenvalues and _get_occupations """
    import re 
    from os.path import exists, join
    from numpy import array
    from ..files import OUTCAR

    path = OUTCAR 
    if len(self.directory): path = join(self.directory, path)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      found = re.compile(r"""k-point\s+(\d+)\s*:\s*(\S+)\s+(\S+)\s+(\S+)$""")
      in_kpoint = -1
      kp_result = []
      for line in file:
        if in_kpoint > 0: 
          data = line.split()
          if len(data) == 3: kp_result.append(float(data[which]))
          else: 
            result.append(kp_result)
            kp_result = []
            in_kpoint = -1
        elif in_kpoint == 0: in_kpoint = 1
        else:
          match = found.search(line)
          if match != None:  
            if int(match.group(1)) == 1: result = []
            in_kpoint = 0
    return array(result, dtype="float64")

  @_bound_mpi_extraction
  def _get_eigenvalues(self):
    """ Returns eigenvalues 

        @return: a two-dimension numpy nxm array of eigenvalues, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(1)
  eigenvalues = property(_get_eigenvalues)
  r""" A matrix of eigenvalues where each row contains eigenvalues of a single kpoint. """

  @_bound_mpi_extraction
  def _get_occupations(self):
    """ Returns band-occupations 

        @return: a two-dimension numpy nxm array of occupations, with n the
                 number of kpoints and m the number of bands.
    """
    return self._get_eigocc(2)
  occupations = property(_get_occupations)
  r""" A matrix of eigenvalues where each row contains band occupations of a single kpoint. """

