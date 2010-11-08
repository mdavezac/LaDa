""" Subpackage containing extraction methods for VASP-DFT data from output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract', 'IOMixin']
from ...opt.decorators import make_cached, broadcast_result
from ...opt import AbstractExtractBase

def dft_method(method):
  """ Marks method as dft method. """
  if hasattr(method, 'fget'): method.fget.__extraction_type__ = 'dft'
  else: method.__extraction_type__ = 'dft'
  return method
def gw_method(method):
  """ Marks method as dft method. """
  if hasattr(method, '__get__'): method.fget.__extraction_type__ = 'gw'
  else: method.__extraction_type__ = 'gw'
  return method

class Extract(AbstractExtractBase):
  """ Implementation class for extracting data from VASP output """

  def __init__(self, directory = None, comm = None):
    """ Initializes the extraction class. 

        :Parameters: 
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory. Can also be the path to the
            OUTCAR file itself. 
          comm : boost.mpi.communicator or None
            MPI group communicator. Extraction will be performed for all procs
            in the group. In serial mode, comm can be None.
    """
    super(Extract, self).__init__(directory, comm)
    
  @dft_method
  @property 
  def contcar_path(self):
    """ Returns path to CONTCAR file.

        :raise IOError: if the CONTCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.CONTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return path
    
  @property
  def exports(self):
    """ List of files to export. """
    from os.path import join, exists
    return [ join(self.directory, u) for u in [self.OUTCAR, self.FUNCCAR] \
             if exists(join(self.directory, u)) ]


  def _search_OUTCAR(self, regex):
    """ Looks for all matches. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    result = []
    regex  = compile(regex)
    with self.__outcar__() as file:
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

    result = []
    regex  = compile(regex)
    with self.__outcar__() as file: lines = file.readlines()
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
  def algo(self):
    """ Returns the kind of algorithms. """
    result = self._find_first_OUTCAR(r"""^\s*ALGO\s*=\s*(\S+)\s*""")
    if result == None: return 'Normal'
    return result.group(1).lower()

  @property
  def is_dft(self):
    """ True if this is a DFT calculation, as opposed to GW. """
    return self.algo not in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
  @property
  def is_gw(self):
    """ True if this is a GW calculation, as opposed to DFT. """
    return self.algo in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
    

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def functional(self):
    """ Returns vasp functional used for calculation.

        Requires file FUNCCAR to be present.
    """
    from os.path import exists, join
    from cPickle import load
    try: path = self.funccar_path
    except IOError: return None
    with self.__funccar__ as file: return load(file)

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
      if not exists(join(self.directory, path)): return False
      
    regex = r"""General\s+timing\s+and\s+accounting\s+informations\s+for\s+this\s+job"""
    return self._find_last_OUTCAR(regex) != None


  @broadcast_result(attr=True, which=0)
  def _starting_structure_data(self):
    """ Structure at start of calculations. """
    from re import compile
    from numpy import array, zeros, dot

    species_in = self.species

    cell = zeros((3,3), dtype="float64")
    atoms = []

    with self.__outcar__() as file: 
      atom_index, cell_index = None, None
      cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
      atom_re = compile(r"""^\s*position\s+of\s+ions\s+in\s+fractional\s+coordinates""")
      for line in file:
        if cell_re.search(line) != None: break
      for i in range(3):
        cell[:,i] = array(file.next().split()[:3], dtype='float64')
      for line in file:
        if atom_re.search(line) != None: break
      for line in file:
        data = line.split()
        if len(data) != 3: break
        atoms.append(dot(cell, array(data, dtype='float64')))

    return cell, atoms

  @property
  @make_cached
  def starting_structure(self):
    """ Structure at start of calculations. """
    from ...crystal import Structure
    from quantities import eV
    cell, atoms = self._starting_structure_data()
    structure = Structure()
    # tries to find adequate name for structure.
    try: name = self.system
    except RuntimeError: name = ''
    structure.name = self.name
    if len(name) == 0 or name == 'POSCAR created by SUPERPOSCAR':
      try: title = self.system
      except RuntimeError: title = ''
      if len(title) != 0: structure.name = title

    structure.energy = 0e0
    structure.cell = cell
    structure.scale = 1e0
    assert len(self.species) == len(self.ions_per_specie),\
           RuntimeError("Number of species and of ions per specie incoherent.")
    assert len(atoms) == sum(self.ions_per_specie),\
           RuntimeError('Number of atoms per specie does not sum to number of atoms.')
    for specie, n in zip(self.species,self.ions_per_specie):
      for i in range(n): structure.add_atom = atoms.pop(0), specie

    if (self.isif == 0 or self.nsw == 0 or self.ibrion == -1) and self.is_dft:
      structure.energy = float(self.total_energy.rescale(eV))
    return structure

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _structure_data(self):
    """ Greps cell and positions from OUTCAR. """
    from re import compile
    from numpy import array, zeros

    species_in = self.species

    cell = zeros((3,3), dtype="float64")
    atoms = []

    with self.__outcar__() as file: lines = file.readlines()

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

    if self.isif == 0 or self.nsw == 0 or self.ibrion == -1:
      return self.starting_structure


    species_in = self.species
    try: cell, atoms = self._structure_data
    except: return self.contcar_structure

    structure = Structure()
    # tries to find adequate name for structure.
    try: name = self.system
    except RuntimeError: name = ''
    structure.name = self.name
    if len(name) == 0 or name == 'POSCAR created by SUPERPOSCAR':
      try: title = self.system
      except RuntimeError: title = ''
      if len(title) != 0: structure.name = title
    
    structure.energy = float(self.total_energy.rescale(eV)) if self.is_dft else 0e0
    structure.cell = array(cell, dtype="float64")
    structure.scale = 1e0
    assert len(self.species) == len(self.ions_per_specie),\
           RuntimeError("Number of species and of ions per specie incoherent.")
    for specie, n in zip(self.species,self.ions_per_specie):
      for i in range(n):
        structure.add_atom = array(atoms.pop(0), dtype="float64"), specie

    return structure

  @property
  @make_cached
  def contcar_structure(self):
    """ Greps structure from CONTCAR. """
    from os.path import exists, join
    from ...crystal import read_poscar
    from quantities import eV

    species_in = self.species

    result = read_poscar(species_in, self.__contcar__(), comm=self.comm)
    structure.energy = float(self.total_energy.rescale(eV)) if self.is_dft else 0e0
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
  def isif(self):
    """ Greps ISIF from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISIF\s*=\s*(-?\d+)\s+""")
    if result == None: return None
    return int(result.group(1))
  
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nsw(self):
    """ Greps NSW from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*NSW\s*=\s*(-?\d+)\s+""")
    if result == None: return None
    return int(result.group(1))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ibrion(self):
    """ Greps IBRION from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*IBRION\s*=\s*(-?\d+)\s+""")
    if result == None: return None
    return int(result.group(1))

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

    result = []
    with self.__outcar__() as file:
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

    result = []
    with self.__outcar__() as file:
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

    with self.__outcar__() as file: lines = file.readlines()
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

    with self.__outcar__() as file: lines = file.readlines()
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

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ionic_charges(self):
    """ Greps ionic_charges from OUTCAR."""
    from numpy import array
    regex = """^\s*ZVAL\s*=\s*(.*)$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find ionic_charges in OUTCAR")
    return array([float(u) for u in result.group(1).split()])

  @property
  @make_cached
  def valence(self):
    """ Greps total energy from OUTCAR."""
    from numpy import array
    ionic = self.ionic_charges
    species = self.species
    atoms = [u.type for u in self.structure.atoms]
    result = 0
    for c, s in zip(ionic, species): result += c * atoms.count(s)
    return result
  
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nelect(self):
    """ Greps nelect from OUTCAR."""
    regex = """^\s*NELECT\s*=\s*(\S+)\s+total\s+number\s+of\s+electrons\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) 

  @property
  @make_cached
  def charge(self):
    """ Greps total charge in the system from OUTCAR."""
    return self.valence-self.nelect

  @dft_method
  @property
  @make_cached
  def charge_corrections(self):
     """ First and Third order charge corrections.
     
         Computes first and third order charge corrections according to `Lany
         and Zunger, PRB 78, 235104 (2008)`__. Calculations are
         done for the correct charge of the system and a static dielectric
         constant epsilon=1. For other static dielectric constants, use:

         >>> correction = output.charge_corrections / epsilon

         For conventional and unit-cells of Ga2MnO4 spinels, the charge
         corrections are converged to roughly 1e-5 eV (for singly charged).

         .. __:  http://dx.doi.org/10.1103/PhysRevB.78.235104
     """
     from ...crystal.defects import charge_corrections
     return charge_corrections( self.structure, charge=self.charge, \
                                epsilon=1e0, n=125, cutoff=15e1 )

  @dft_method
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

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def energies_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(2)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def all_total_energies(self):
    """ Greps total energies for all electronic steps from OUTCAR."""
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy =\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @dft_method
  @property
  def cbm(self):
    """ Returns Condunction Band Minimum. """
    from numpy import min
    if self.ispin == 2:
      assert 2 * self.eigenvalues.shape[2] > self.valence + 2,\
             RuntimeError("Not enough bands were computed.")
      return min(self.eigenvalues[:, :, self.valence/2])
    else:
      assert 2 * self.eigenvalues.shape[1] > (self.valence/2) + 1,\
             RuntimeError("Not enough bands were computed.")
      return min(self.eigenvalues[:, self.valence/2])

  @dft_method
  @property
  def vbm(self):
    """ Returns Valence Band Maximum. """
    from numpy import max
    if self.ispin == 2:
      assert 2 * self.eigenvalues.shape[2] > self.valence,\
             RuntimeError("Not enough bands were computed.")
      return max(self.eigenvalues[:, :, self.valence/2-1])
    else:
      assert 2 * self.eigenvalues.shape[1] > self.valence,\
             RuntimeError("Not enough bands were computed.")
      return max(self.eigenvalues[:, self.valence/2-1])

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def total_energies(self):
    """ Greps total energies for all ionic steps from OUTCAR."""
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def total_energy(self):
    """ Greps total energy from OUTCAR."""
    from quantities import eV
    regex = """energy\s+without\s+entropy=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) * eV

  @dft_method
  @property
  def energy(self): 
    """ Alias for total_energy. """
    return self.total_energy

  @dft_method
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

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def moment(self):
    """ Returns magnetic moment from OUTCAR. """
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find magnetic moment in OUTCAR")
    return float(result.group(2))

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nb_electrons(self):
    """ Returns magnetic moment from OUTCAR. """
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find number of electrons in OUTCAR")
    return float(result.group(1))

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def pressures(self):
    """ Greps all pressures from OUTCAR """
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find pressures in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find pressures in OUTCAR")
    return result * kB

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def pressure(self):
    """ Greps last pressure from OUTCAR """
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find pressure in OUTCAR")
    return float(result.group(1)) * kB

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def alphabet(self):
    """ Greps alpha+bet from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find alpha+bet in OUTCAR")
    return float(result.group(3))

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def xc_g0(self):
    """ Greps XC(G=0) from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find xc(G=0) in OUTCAR")
    return float(result.group(2))

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def pulay_pressure(self):
    from quantities import kbar as kB
    """ Greps pressure from OUTCAR """
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find pulay pressure in OUTCAR")
    return float(result.group(2)) * kB

  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    from os.path import exists, join
    from re import compile, search, X as re_X

    result = None
    with self.__outcar__() as file:

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

    result = []
    with self.__outcar__() as file: lines = file.readlines()
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

  @dft_method
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

  @dft_method
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


  @dft_method
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

  @dft_method
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



  @dft_method
  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def electropot(self):
    """ Greps average atomic electrostatic potentials from OUTCAR. """
    from os.path import exists, join
    from re import compile, X as reX
    from numpy import array
    from quantities import eV

    with self.__outcar__() as file: lines = file.readlines()
    regex = compile(r"""average\s+\(electrostatic\)\s+potential\s+at\s+core""", reX)
    for i, line in enumerate(lines[::-1]):
      if regex.search(line) != None: break
    assert -i + 2 < len(lines), RuntimeError("Could not find average atomic potential in file.")
    result = []
    for line in lines[-i+2:]:
      data = line.split()
      if len(data) == 0: break
      result.extend( [float(u) for i, u in enumerate(data) if i % 2 == 1] )
        
    return array(result, dtype="float64") * eV


  @gw_method
  @property
  @make_cached
  def dft_eigenvalues(self):
    """ Greps DFT eigenvalues from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(1), dtype="float64") * eV

  @gw_method
  @property
  @make_cached
  def qp_eigenvalues(self):
    """ Greps quasi-particle eigenvalues from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(2), dtype="float64") * eV

  @gw_method
  @property
  @make_cached
  def self_energies(self):
    """ Greps self-energies of each eigenvalue from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(3), dtype="float64") * eV

  @gw_method
  @property
  @make_cached
  def occupations(self):
    """ Greps occupations from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    from quantities import eV
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(3), dtype="float64") * eV


class IOMixin(object):
  """ A mixin base clase which controls file IO. 

      Defines special property with file-like behaviors. 
      Makes it easier to change the behavior of the extraction class.
  """
  def __init__(self, directory=None, OUTCAR=None, FUNCCAR=None, CONTCAR=None):
    """ Initializes the extraction class. 

        :Parameters: 
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory. Can also be the path to the
            OUTCAR file itself. 
          OUTCAR : str or None
            If given, this name will be used, rather than files.OUTCAR.
          CONTCAR : str or None
            If given, this name will be used, rather than files.CONTCAR.
          FUNCCAR : str or None
            If given, this name will be used, rather than files.FUNCCAR.
    """
    from .. import files
    
    super(IOMixin, self).__init__()

    self.OUTCAR  = OUTCAR if OUTCAR != None else files.OUTCAR
    """ Filename of the OUTCAR file from VASP. """
    self.CONTCAR  = CONTCAR if CONTCAR != None else files.CONTCAR
    """ Filename of the CONTCAR file from VASP. """
    self.FUNCCAR  = FUNCCAR if FUNCCAR != None else files.FUNCCAR
    """ Filename of the FUNCCAR file containing the pickled functional. """

  def __outcar__(self):
    """ Returns path to OUTCAR file.

        :raise IOError: if the OUTCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  def __funccar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.FUNCCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  def __contcar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.CONTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')
