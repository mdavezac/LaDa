""" Subpackage containing extraction methods for VASP common to both GW and DFT. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from quantities import g, cm, eV
from ...opt.decorators import make_cached, broadcast_result
from ...opt.json import array as json_array, structure as json_structure, section as json_section, \
                        unit as json_unit

class Extract(object):
  """ Implementation class for extracting data from VASP output """

  def __init__(self):
    """ Initializes the extraction class. """
    object.__init__(self)

  @property 
  @make_cached
  @broadcast_result(attr=True, which=0)
  def algo(self):
    """ Returns the kind of algorithms. """
    result = self._find_first_OUTCAR(r"""^\s*ALGO\s*=\s*(\S+)\s*""")
    if result is None: return 'Normal'
    return result.group(1).lower()

  @property
  @json_section("input")
  def is_dft(self):
    """ True if this is a DFT calculation, as opposed to GW. """
    try: return self.algo not in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
    except: return False
  @property
  @json_section("input")
  def is_gw(self):
    """ True if this is a GW calculation, as opposed to DFT. """
    try: return self.algo in ['gw', 'gw0', 'chi', 'scgw', 'scgw0'] 
    except: return False
    
  @property 
  @json_section("input")
  @json_unit(eV)
  @broadcast_result(attr=True, which=0)
  def encut(self):
    """ Energy cutoff. """
    return float(self._find_first_OUTCAR(r"ENCUT\s*=\s*(\S+)").group(1)) * eV


  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def functional(self):
    """ Returns vasp functional used for calculation.

        Requires file FUNCCAR to be present.
    """
    from cPickle import load
    with self.__funccar__() as file: return load(file)

  @property
  def vasp(self):
    """ Deprecated. """
    from warnings import warn
    warn(DeprecationWarning('vasp attribute is deprecated in favor of functional.'), stacklevel=2)
    return self.functional

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def success(self):
    """ Checks that VASP run has completed. 

        At this point, checks for the existence of OUTCAR.
        Then checks that timing stuff is present at end of OUTCAR.
    """
    from os.path import exists, join

    if hasattr(self, "OUTCAR"): 
      for path in [self.OUTCAR]:
        if not exists(join(self.directory, path)): return False
      
    regex = r"""General\s+timing\s+and\s+accounting\s+informations\s+for\s+this\s+job"""
    return self._find_last_OUTCAR(regex) is not None


  @broadcast_result(attr=True, which=0)
  def _starting_structure_data(self):
    """ Structure at start of calculations. """
    from re import compile
    from numpy import array, zeros, dot
    from numpy.linalg import inv

    cell = zeros((3,3), dtype="float64")
    atoms = []

    with self.__outcar__() as file: 
      atom_index, cell_index = None, None
      cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
      atom_re = compile(r"""^\s*position\s+of\s+ions\s+in\s+fractional\s+coordinates""")
      for line in file:
        if cell_re.search(line) is not None: break
      data = []
      for i in range(3):
        data.append(file.next().split())
      try: 
        for i in range(3): cell[:,i] = array(data[i][:3], dtype='float64')
      except: 
        for i in range(3): cell[i, :] = array(data[i][-3:], dtype='float64')
        cell = inv(cell)
      for line in file:
        if atom_re.search(line) is not None: break
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
    assert len(self.species) == len(self.stoechiometry),\
           RuntimeError("Number of species and of ions per specie incoherent.")
    assert len(atoms) == sum(self.stoechiometry),\
           RuntimeError('Number of atoms per specie does not sum to number of atoms.')
    for specie, n in zip(self.species,self.stoechiometry):
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
    from numpy.linalg import inv
    from numpy import array, zeros

    cell = zeros((3,3), dtype="float64")
    atoms = []

    with self.__outcar__() as file: lines = file.readlines()

    atom_index, cell_index = None, None
    atom_re = compile(r"""^\s*POSITION\s+""")
    cell_re = compile(r"""^\s*direct\s+lattice\s+vectors\s+""")
    for index, line in enumerate(lines[::-1]):
      if atom_re.search(line) is not None: atom_index = index - 1
      if cell_re.search(line) is not None: cell_index = index; break
    assert atom_index is not None and cell_index is not None,\
           RuntimeError("Could not find structure description in OUTCAR.")
    try: 
      for i in range(3): cell[:,i] = array(lines[-cell_index+i].split()[:3], dtype="float64")
    except: 
      for i in range(3): cell[i,:] = array(lines[-cell_index+i].split()[-3:], dtype="float64")
      cell = inv(cell)
    for i in range(3):
      cell[:,i] = [float(u) for u in lines[-cell_index+i].split()[:3]]
    while atom_index > 0 and len(lines[-atom_index].split()) == 6:
      atoms.append( array([float(u) for u in lines[-atom_index].split()[:3]], dtype="float64") )
      atom_index -= 1

    return cell, atoms

  @property
  @json_section(None)
  @json_structure
  @make_cached
  def structure(self):
    """ Greps structure and total energy from OUTCAR. """
    from numpy import array
    from quantities import eV
    from ...crystal import Structure

    if self.isif == 0 or self.nsw == 0 or self.ibrion == -1:
      return self.starting_structure


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
    assert len(self.species) == len(self.stoechiometry),\
           RuntimeError("Number of species and of ions per specie incoherent.")
    for specie, n in zip(self.species,self.stoechiometry):
      for i in range(n):
        structure.add_atom = array(atoms.pop(0), dtype="float64"), specie

    return structure

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def LDAUType(self):
    """ Type of LDA+U performed. """
    type = self._find_first_OUTCAR(r"""LDAUTYPE\s*=\s*(\d+)""")
    if type == None: return 0
    type = int(type.group(1))
    if type == 1: return "liechtenstein"
    elif type == 2: return "dudarev"
    return type

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def HubbardU_NLEP(self):
    """ Hubbard U/NLEP parameters. """
    from ..specie import U as ldaU, nlep
    from re import M
    type = self._find_first_OUTCAR(r"""LDAUTYPE\s*=\s*(\d+)""")
    if type == None: return {}
    type = int(type.group(1))

    species = tuple([ u.group(1) for u in self._search_OUTCAR(r"""VRHFIN\s*=\s*(\S+)\s*:""") ])

    # first look for standard VASP parameters.
    groups = r"""\s*((?:(?:\+|-)?\d+(?:\.\d+)?\s*)+)\s*\n"""
    regex = r"""\s*angular\s+momentum\s+for\s+each\s+species\s+LDAUL\s+=""" + groups \
          + r"""\s*U\s+\(eV\)\s+for\s+each\s+species\s+LDAUU\s+="""         + groups \
          + r"""\s*J\s+\(eV\)\s+for\s+each\s+species\s+LDAUJ\s+="""         + groups
    result = {k: [] for k in species}
    for found in self._search_OUTCAR(regex, M):
      moment = found.group(1).split()
      LDAU   = found.group(2).split()
      LDAJ   = found.group(3).split()
      for specie, m, U, J in zip(species, moment, LDAU, LDAJ):
        if int(m) != -1: result[specie].append(ldaU(type, int(m), float(U), float(J)))
    for key in result.keys():
      if len(result[key]) == 0: del result[key]
    if len(result) != 0: return result

    # then look for nlep parameters.
    regex = r"""\s*angular\s+momentum\s+for\s+each\s+species,\s+LDAU#\s*(?:\d+)\s*:\s*L\s*=""" + groups \
          + r"""\s*U\s+\(eV\)\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*U\s*=""" + groups \
          + r"""\s*J\s+\(eV\)\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*J\s*=""" + groups \
          + r"""\s*nlep\s+for\s+each\s+species,\s+LDAU\#\s*(?:\d+)\s*:\s*O\s*=""" + groups 
    result = {k: [] for k in species}
    for found in self._search_OUTCAR(regex, M):
      moment = found.group(1).split()
      LDAU   = found.group(2).split()
      LDAJ   = found.group(3).split()
      funcs  = found.group(4).split()
      for specie, m, U, J, func in zip(species, moment, LDAU, LDAJ, funcs):
        if int(m) == -1: continue
        if int(func) == 1: result[specie].append(ldaU(type, int(m), float(U), float(J)))
        else: result[specie].append(nlep(type, int(m), float(U), float(J)))
    for key in result.keys():
      if len(result[key]) == 0: del result[key]
    return result

  @property
  @json_section("input")
  @make_cached
  def pseudopotential(self):
    """ Title of the first POTCAR. """
    return self._find_first_OUTCAR(r"""POTCAR:.*""").group(0).split()[1]


  @property
  @json_section("input")
  @make_cached
  def volume(self): 
    """ Unit-cell volume. """
    from numpy import abs
    from numpy.linalg import det
    from quantities import angstrom
    return abs(self.structure.scale * det(self.structure.cell)) * angstrom**3

  @property 
  @json_section("input")
  @make_cached
  def reciprocal_volume(self):
    """ Reciprocal space volume (including 2pi). """
    from numpy import abs, pi
    from numpy.linalg import det, inv
    from quantities import angstrom
    return abs(det(inv(self.structure.scale * self.structure.cell))) * (2e0*pi/angstrom)**3

  @property
  @json_section("output")
  @json_unit(g/cm**3)
  @make_cached
  def density(self):
    """ Computes density of the material. """
    from quantities import atomic_mass_unit as amu
    from ... import periodic_table as pt
    result = 0e0 * amu;
    for atom, n in zip(self.species, self.stoechiometry): 
      if atom not in pt.__dict__: return None;
      result += pt.__dict__[atom].atomic_weight * float(n) * amu 
    return (result / self.volume).rescale(g/cm**3)

  @property
  @make_cached
  def contcar_structure(self):
    """ Greps structure from CONTCAR. """
    from ...crystal import read_poscar
    from quantities import eV

    species_in = self.species

    result = read_poscar(species_in, self.__contcar__(), comm=self.comm)
    result.energy = float(self.total_energy.rescale(eV)) if self.is_dft else 0e0
    return result

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def stoechiometry(self):
    """ Stoechiometry of the compound. """
    result = self._find_first_OUTCAR(r"""\s*ions\s+per\s+type\s*=.*$""")
    if result is None: return None
    return [int(u) for u in result.group(0).split()[4:]]
  def ions_per_specie(self): 
    """ Alias for stoechiometry. """
    return self.stoechiometry()

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def species(self):
    """ Greps species from OUTCAR. """
    return tuple([ u.group(1) for u in self._search_OUTCAR(r"""VRHFIN\s*=\s*(\S+)\s*:""") ])

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def isif(self):
    """ Greps ISIF from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*ISIF\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))
  
  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nsw(self):
    """ Greps NSW from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*NSW\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ibrion(self):
    """ Greps IBRION from OUTCAR. """
    result = self._find_first_OUTCAR(r"""\s*IBRION\s*=\s*(-?\d+)\s+""")
    if result is None: return None
    return int(result.group(1))

  @property
  @json_section("input")
  @json_array("float64")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def kpoints(self):
    """ Greps k-points from OUTCAR.
    
        Numpy array where each row is a k-vector in cartesian units. 
    """
    from re import compile
    from numpy import array

    result = []
    found_generated = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
    found_read = compile(r"""k-points in units of 2pi/SCALE and weight""")
    with self.__outcar__() as file:
      found = 0
      for line in file:
        if found_generated.search(line) is not None: found = 1; break
        elif found_read.search(line) is not None: found = 2; break
      if found == 1:
        found = compile(r"""Following\s+cartesian\s+coordinates:""")
        for line in file:
          if found.search(line) is not None: break
        file.next()
        for line in file:
          data = line.split()
          if len(data) != 4: break;
          result.append( data[:3] )
        return array(result, dtype="float64") 
      if found == 2:
        for line in file:
          data = line.split()
          if len(data) == 0: break
          result.append(data[:3])
        return array(result, dtype="float64") 

  @property
  @json_section("input")
  @json_array("float64")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def multiplicity(self):
    """ Greps multiplicity of each k-point from OUTCAR. """
    from re import compile
    from numpy import array

    result = []
    # case where kpoints where generated by vasp.
    found_generated = compile(r"""Found\s+(\d+)\s+irreducible\s+k-points""")
    # case where kpoints where given by hand.
    found_read = compile(r"""k-points in units of 2pi/SCALE and weight""")
    with self.__outcar__() as file:
      found = 0
      for line in file:
        if found_generated.search(line) is not None: found = 1; break
        elif found_read.search(line) is not None: found = 2; break
      if found == 1:
        found = compile(r"""Following\s+cartesian\s+coordinates:""")
        for line in file:
          if found.search(line) is not None: break
        file.next()
        for line in file:
          data = line.split()
          if len(data) != 4: break;
          result.append( data[-1] )
        return array(result, dtype="float64") 
      elif found == 2:
        for line in file:
          data = line.split()
          if len(data) == 0: break
          result.append(float(data[3]))
        return array([round(r*float(len(result))) for r in result], dtype="float64")

  @property 
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ispin(self):
    """ Greps ISPIN from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*ISPIN\s*=\s*(1|2)\s+""")
    assert result is not None, RuntimeError("Could not extract ISPIN from OUTCAR.")
    return int(result.group(1))

  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def name(self):
    """ Greps POSCAR title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*POSCAR\s*=.*$""")
    assert result is not None, RuntimeError("Could not extract POSCAR title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:]
    return result.rstrip().lstrip()

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def system(self):
    """ Greps system title from OUTCAR. """
    result = self._find_first_OUTCAR(r"""^\s*SYSTEM\s*=.*$""")
    assert result is not None, RuntimeError("Could not extract SYSTEM title from OUTCAR.")
    result = result.group(0)
    result = result[result.index('=')+1:].rstrip().lstrip()
    if result[0] == '"': result = result[1:]
    if result[-1] == '"': result = result[:-1]
    return result

  @broadcast_result(attr=True, which=0)
  def _unpolarized_values(self, which):
    """ Returns spin-unpolarized eigenvalues and occupations. """
    from re import compile, finditer
    import re

    with self.__outcar__() as file: lines = file.readlines()
    # Finds last first kpoint.
    spin_comp1_re = compile(r"\s*k-point\s+1\s*:\s*(\S+)\s+(\S+)\s+(\S+)\s*")
    found = None
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found is not None: break
    assert found is not None, RuntimeError("Could not extract eigenvalues/occupation from OUTCAR.")

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
    for kp in finditer(kp_re, "".join(lines[-i-1:]), re.M):
      dummy = [u.split() for u in kp.group(0).split('\n')[skip:]]
      results.append([float(u[which]) for u in dummy if len(u) == cols])
    return results

  @broadcast_result(attr=True, which=0)
  def _spin_polarized_values(self, which):
    """ Returns spin-polarized eigenvalues and occupations. """
    from re import compile, finditer
    import re

    with self.__outcar__() as file: lines = file.readlines()
    # Finds last spin components.
    spin_comp1_re = compile(r"""\s*spin\s+component\s+(1|2)\s*$""")
    spins = [None,None]
    for i, line in enumerate(lines[::-1]):
      found = spin_comp1_re.match(line)
      if found is None: continue
      if found.group(1) == '1': 
        assert spins[1] is not None, \
               RuntimeError("Could not find two spin components in OUTCAR.")
        spins[0] = i
        break
      else:  spins[1] = i
    assert spins[0] is not None and spins[1] is not None,\
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
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def ionic_charges(self):
    """ Greps ionic_charges from OUTCAR."""
    regex = """^\s*ZVAL\s*=\s*(.*)$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find ionic_charges in OUTCAR")
    return [float(u) for u in result.group(1).split()]

  @property
  @json_section("input")
  @make_cached
  def valence(self):
    """ Greps number of valence bands from OUTCAR."""
    ionic = self.ionic_charges
    species = self.species
    atoms = [u.type for u in self.structure.atoms]
    result = 0
    for c, s in zip(ionic, species): result += c * atoms.count(s)
    return result
  
  @property
  @json_section("input")
  @make_cached
  @broadcast_result(attr=True, which=0)
  def nelect(self):
    """ Greps nelect from OUTCAR."""
    regex = """^\s*NELECT\s*=\s*(\S+)\s+total\s+number\s+of\s+electrons\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) 

  @property
  @json_section("input")
  @make_cached
  def charge(self):
    """ Greps total charge in the system from OUTCAR."""
    return self.valence-self.nelect

  @property
  @make_cached
  def nbands(self):
    """ Number of bands in calculation. """
    result = self._find_first_OUTCAR("""NBANDS\s*=\s*(\d+)""")
    assert result is not None, RuntimeError("Could not find NBANDS in OUTCAR.")
    return int(result.group(1))


  def iterfiles(self, **kwargs):
    """ iterates over input/output files. 
    
        :kwarg errors: Include stderr files.
        :type errors: bool
        :kwarg incar: Include INCAR file
        :type incar: bool
        :kwarg wavecar: Include WAVECAR file
        :type wavecar: bool
        :kwarg doscar: Include CHGCAR file
        :type doscar: bool
        :kwarg chgcar: Include CHGCAR file
        :type chgcar: bool
        :kwarg poscar: Include POSCAR file
        :type poscar: bool
        :kwarg contcar: Include CONTCAR file
        :type contcar: bool
        :kwarg procar: Include PROCAR file
        :type procar: bool
    """
    from os.path import exists, join
    from glob import iglob
    from itertools import chain
    files = [self.OUTCAR, self.FUNCCAR]
    try: files.append(self.functional.STDOUT)
    except: pass
    if kwargs.get('errors', False): 
      try: files.append(self.functional.STDERR)
      except: pass
    if kwargs.get('incar', False):   files.append('INCAR')
    if kwargs.get('wavecar', False): files.append('WAVECAR')
    if kwargs.get('doscar', False):  files.append('DOSCAR')
    if kwargs.get('chgcar', False):  files.append('CHGCAR')
    if kwargs.get('poscar', False):  files.append('POSCAR')
    if kwargs.get('contcar', False): files.append('CONTCAR')
    if kwargs.get('procar', False): files.append('PROCAR')
    for file in files:
      file = join(self.directory, file)
      if exists(file): yield file
    # Add RelaxCellShape directories.
    for dir in chain( iglob(join(self.directory, "relax_cellshape/[0-9]/")),
                      iglob(join(self.directory, "relax_cellshape/[0-9][0-9]/")),
                      iglob(join(self.directory, "relax_ions/[0-9]/")),
                      iglob(join(self.directory, "relax_ions/[0-9][0-9]/")) ):
      a = self.__class__(dir, comm=self.comm)
      for file in a.iterfiles(**kwargs): yield file
