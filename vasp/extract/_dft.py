""" Subpackage containing extraction methods for VASP-DFT data from output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ...opt.decorators import make_cached, broadcast_result
class _ExtractImpl(object):
  """ Implementation class for extracting data from VASP output """

  def __init__(self, directory = None, comm = None, OUTCAR = None):
    """ Initializes the extraction class. 

        :Parameters: boost.mpi.Communicator
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory.
          comm : boost.mpi.communicator or None
            MPI group communicator. Extraction will be performed for all procs
            in the group. In serial mode, comm can be None.
          OUTCAR : str or None
            Name of the OUTCAR file.
    """
    from os import getcwd
    from .. import files
    from ...opt import RelativeDirectory
    super(_ExtractImpl, self).__init__()

    if directory == None: directory = getcwd()
    self._directory = RelativeDirectory(directory, hook=self.uncache)
    """ Directory where to check for output. """
    self.comm = comm
    """ MPI group communicator. """
    self.OUTCAR = files.OUTCAR if OUTCAR == None else OUTCAR
    """ Filename of the OUTCAR file from VASP. """
    self.CONTCAR = files.CONTCAR
    """ Filename of the CONTCAR file from VASP. """
    self.FUNCCAR = files.FUNCCAR
    """ Filename of the FUNCCAR file containing the pickled functional. """
    
  @property
  def exports(self):
    """ List of files to export. """
    from os.path import join, exists
    return [ join(self.directory, u) for u in [self.OUTCAR, self.FUNCCAR] \
             if exists(join(self.directory, u)) ]

  @property
  def directory(self):
    """ Directory with VASP output files """
    return self._directory.path
  @directory.setter
  def directory(self, value): self._directory.path = value

  def uncache(self): 
    """ Removes cached results.

        After this outputs are re-read from file.
    """
    from ...opt.decorators import uncache
    uncache(self)

  def solo(self):
    """ Extraction on a single process.

        Sometimes, it is practical to perform extractions on a single process
        only, eg without blocking mpi calls. ``self.solo()`` returns an
        extractor for a single process:
        
        >>> # prints only on proc 0.
        >>> if boost.mpi.world.rank == 0: print extract.solo().structure
    """
    from copy import deepcopy
    return deepcopy(self) if self.comm != None else self

  def __getstate__(self):
    from os.path import relpath
    d = self.__dict__.copy()
    if "_directory" in d: d["_directory"].hook = None
    if "comm" in d: del d["comm"]
    return d
  def __setstate__(self, arg):
    self.__dict__.update(arg)
    self.comm = None
    if hasattr(self, "_directory"): self._directory.hook = self.uncache


  def _search_OUTCAR(self, regex):
    """ Looks for all matches. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    regex  = compile(regex)
    with open(path, "r") as file:
      for line in file: 
        found = regex.search(line)
        if found != None: yield found

  def _find_first_OUTCAR(self, regex):
    """ Returns first result from a regex. """
    for first in self._search_OUTCAR(regex): return first
    return None

  def _rsearch_OUTCAR(self, regex):
    """ Looks for all matches starting from the end. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    regex  = compile(regex)
    with open(path, "r") as file: lines = file.readlines()
    for line in lines[::-1]:
      found = regex.search(line)
      if found != None: yield found

  def _find_last_OUTCAR(self, regex):
    """ Returns first result from a regex. """
    for last in self._rsearch_OUTCAR(regex): return last
    return None

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def functional(self):
    """ Returns vasp functional used for calculation.

        Requires file L{FUNCCAR} to be present.
    """
    from os.path import exists, join
    from cPickle import load
    path = self.FUNCCAR if len(self.directory) == 0 else join(self.directory, self.FUNCCAR)
    if not exists(path): return None
    with open(path, "r") as file: return load(file)

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks that VASP run has completed. 

        At this point, checks for the existence of OUTCAR.
        Then checks that timing stuff is present at end of OUTCAR.
    """
    from os.path import exists, join
    import re

    for path in [self.OUTCAR]:
      if self.directory != "": path = join(self.directory, path)
      if not exists(path): return False
      
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)

    regex = r"""General\s+timing\s+and\s+accounting\s+informations\s+for\s+this\s+job"""
    return self._find_last_OUTCAR(regex) != None

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _structure_data(self):
    """ Greps cell and positions from OUTCAR. """
    from os.path import exists, join
    from re import compile
    from numpy import array, zeros

    species_in = self.species

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    cell = zeros((3,3), dtype="float64")
    atoms = []

    with open(path, "r") as file: lines = file.readlines()

    atom_index, cell_index = None, None
    atom_re = compile(r"""^\s*POSITION\s+""")
    cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
    for index, line in enumerate(lines[::-1]):
      if atom_re.search(line) != None: atom_index = index - 1
      if cell_re.search(line) != None: cell_index = index; break
    assert atom_index != None and cell_index != None,\
           RuntimeError("Could not find structure description in OUTCAR.")
    for i in range(3):
      cell[:,i] = [float(u) for u in lines[-cell_index+i].split()[:3]]
    while atom_index > 0 and len(lines[-atom_index].split()) == 6:
      atoms.append( array([float(u) for u in lines[-atom_index].split()[:3]], dtype="float64") )
      atom_index -= 1

    return cell, atoms

  @property
  @make_cached
  def structure(self):
    """ Greps structure and total energy from OUTCAR. """
    from re import compile
    from numpy import array, zeros
    from quantities import eV
    from ...crystal import Structure

    species_in = self.species
    try: cell, atoms = self._structure_data
    except: return self.contcar_structure

    structure = Structure()
    structure.name = self.name
    structure.energy = float(self.total_energy.rescale(eV))
    structure.cell = cell
    structure.scale = 1e0
    assert len(self.species) == len(self.ions_per_specie),\
           RuntimeError("Number of species and of ions per specie incoherent.")
    for specie, n in zip(self.species[::-1],self.ions_per_specie[::-1]):
      for i in range(n):
        structure.add_atom = atoms.pop(), specie

    return structure

  @property
  @make_cached
  def contcar_structure(self):
    """ Greps structure from CONTCAR. """
    from os.path import exists, join
    from ...crystal import read_poscar
    from quantities import eV

    species_in = self.species

    path = self.CONTCAR if len(self.directory) == 0 else join(self.directory, self.CONTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)
    result = read_poscar(species_in, path, comm=self.comm)
    result.energy = float(self.energy.rescale(eV))
    return result

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ions_per_specie(self):
    """ Greps species from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ions\s+per\s+type\s*=.*$""")
    if result == None: return None
    return [int(u) for u in result.group(0).split()[4:]]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def species(self):
    """ Greps species from OUTCAR. """
    return tuple([ u.group(1) for u in self._search_OUTCAR(r"""VRHFIN\s*=\s*(\S+)\s*:""") ])


  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def kpoints(self):
    """ Greps k-points from OUTCAR.
    
        Numpy array where each row is a k-vector in cartesian units. 
    """
    from os.path import exists, join
    from re import compile, search 
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
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

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def multiplicity(self):
    """ Greps multiplicity of each k-point from OUTCAR. """
    from os.path import exists, join
    from re import compile, search 
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
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

  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ispin(self):
    """ Greps ISPIN from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*ISPIN\s*=\s*(1|2)\s+""")
    assert result != None, RuntimeError("Could not extract ISPIN from OUTCAR.")
    return int(result.group(1))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def name(self):
    """ Greps POSCAR title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*POSCAR\s*=.*$""")
    assert result != None, RuntimeError("Could not extract POSCAR title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:]
    return result.rstrip().lstrip()

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def system(self):
    """ Greps system title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*SYSTEM\s*=.*$""")
    assert result != None, RuntimeError("Could not extract SYSTEM title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:]
    return result.rstrip().lstrip()

  @broadcast_result(attr=True, which=0)
  def _unpolarized_values(self, which):
    """ Returns spin-unpolarized eigenvalues and occupations. """
    from re import compile, search, finditer
    import re
    from os.path import exists, join
    from numpy import array

    # checks for existence.
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    with open(path, "r") as file: lines = file.readlines()
    # Finds last first kpoint.
    spin_comp1_re = compile(r"\s*k-point\s+1\s*:\s*(\S+)\s+(\S+)\s+(\S+)\s*")
    found = None
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found != None: break
    assert found != None, RuntimeError("Could not extract eigenvalues/occupation from OUTCAR.")

    # now greps actual results.
    if self.is_dft:
      kp_re = r"\s*k-point\s+(?:\d+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+band\s+energies\s+occupation\s*\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 2, 3
    else: 
      kp_re = r"\s*k-point\s+(?:\d+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+.*\n\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)"\
              r"\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 3, 8
    results = []
    for kp in finditer(kp_re, "".join(lines[-i:]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results.append([float(u[which]) for u in dummy if len(u) == cols])
    return results

  @broadcast_result(attr=True, which=0)
  def _spin_polarized_values(self, which):
    """ Returns spin-polarized eigenvalues and occupations. """
    from re import compile, search, finditer
    import re
    from os.path import exists, join
    from numpy import array

    # checks for existence.
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    with open(path, "r") as file: lines = file.readlines()
    # Finds last spin components.
    spin_comp1_re = compile(r"""\s*spin\s+component\s+(1|2)\s*$""")
    spins = [None,None]
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found == None: continue
      if found.group(1) == '1': 
        assert spins[1] != None, \
               RuntimeError("Could not find two spin components in OUTCAR.")
        spins[0] = i
        break
      else:  spins[1] = i
    assert spins[0] != None and spins[1] != None,\
           RuntimeError("Could not extract eigenvalues/occupation from OUTCAR.")

    # now greps actual results.
    if self.is_dft:
      kp_re = r"\s*k-point\s+(?:\d+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+band\s+energies\s+occupation\s*\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 2, 3
    else: 
      kp_re = r"\s*k-point\s+(?:\d+)\s*:\s*(?:\S+)\s*(?:\S+)\s*(?:\S+)\n"\
              r"\s*band\s+No\.\s+.*\n\n"\
              r"(\s*(?:\d+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)"\
              r"\s+(?:\S+)\s+(?:\S+)\s+(?:\S+)\s*\n)+"
      skip, cols = 3, 8
    results = [ [], [] ]
    for kp in finditer(kp_re, "".join(lines[-spins[0]:-spins[1]]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results[0].append([float(u[which]) for u in dummy if len(u) == cols])
    for kp in finditer(kp_re, "".join(lines[-spins[1]:]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results[1].append([u[which] for u in dummy if len(u) == cols])
    return results

class Extract(_ExtractImpl): 
  """ Extracts output from OUTCAR, including DFT specific stuff. """

  def __init__(self, directory = None, comm = None, OUTCAR = None): 
    """ Initializes the extraction class. 

        :Parameters: boost.mpi.Communicator
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory.
          comm : boost.mpi.communicator or None
            MPI group communicator. Extraction will be performed for all procs
            in the group. In serial mode, comm can be None.
            OUTCAR : str or None
            Name of the OUTCAR file.
    """
    super(Extract, self).__init__(directory, comm, OUTCAR)

  @property
  def is_dft(self): return True

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def energy_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find sigma0 energy in OUTCAR")
    return float(result.group(2)) * eV

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def total_energy(self):
    """ Greps total energy from OUTCAR."""
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) * eV

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def free_energy(self):
    """ Greps total free energy from OUTCAR. """
    from quantities import eV
    regex = r"""free\s+energy\s+TOTEN\s*=\s*(\S+)\s+eV""" 
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find free energy in OUTCAR")
    return float(result.group(1)) * eV

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def fermi_energy(self):
    """ Greps fermi energy from OUTCAR. """
    from quantities import eV
    regex = r"""E-fermi\s*:\s*(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find fermi energy in OUTCAR")
    return float(result.group(1)) * eV

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def moment(self):
    """ Returns magnetic moment from OUTCAR. """
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find magnetic moment in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nb_electrons(self):
    """ Returns magnetic moment from OUTCAR. """
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find number of electrons in OUTCAR")
    return float(result.group(1))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def pressure(self):
    """ Greps pressure from OUTCAR """
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find pressure in OUTCAR")
    return float(result.group(1))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def pulay_pressure(self):
    """ Greps pressure from OUTCAR """
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find pulay pressure in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    from os.path import exists, join
    from re import compile, search, X as re_X

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
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

  def _get_partial_charges_magnetization(self, grep):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included. Implementation also used for magnetization.
    """
    import re 
    from os.path import exists, join
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file: lines = file.readlines()
    found = re.compile(grep) 
    for index in xrange(1, len(lines)+1):
      if found.search(lines[-index]) != None: break 
    if index == len(lines): return None
    index -= 4
    line_re = re.compile(r"""^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$""")
    for i in xrange(0, index): 
      match = line_re.match(lines[-index+i])
      if match == None: break
      result.append([float(match.group(j)) for j in range(1, 5)])
    return array(result, dtype="float64")

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def partial_charges(self):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    return self._get_partial_charges_magnetization(r"""\s*total\s+charge\s*$""")

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def magnetization(self):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    return self._get_partial_charges_magnetization(r"""^\s*magnetization\s*\(x\)\s*$""")


  @property
  @make_cached
  def eigenvalues(self):
    """ Greps eigenvalues from OUTCAR. 

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(1), dtype="float64") * eV

  @property
  @make_cached
  def occupations(self):
    """ Greps occupations from OUTCAR. 

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(2), dtype="float64")
    return array(self._unpolarized_values(2), dtype="float64") 



  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def electropot(self):
    """ Greps average atomic electrostatic potentials from OUTCAR. """
    from os.path import exists, join
    import re
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError("File %s does not exist.\n" % (path))

    with open(path, "r") as file: lines = file.readlines()
    regex = re.compile(r"""average\s+\(electrostatic\)\s+potential\s+at\s+core""", re.X)
    for i, line in enumerate(lines[::-1]):
      if regex.search(line) != None: break
    assert i + 3 < len(lines), RuntimeError("Could not find average atomic potential in file.")
    result = []
    for line in lines[-i+3:]:
      data = line.split()
      if len(data) == 0: break
      result.extend( [float(u) for i, u in enumerate(data) if i % 2 == 1] )
        
    return array(result, dtype="float64")
