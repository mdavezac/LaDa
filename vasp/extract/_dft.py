""" Subpackage containing extraction methods for VASP-DFT data from output. """
from ...opt.decorators import make_cached, broadcast_result
class _ExtractImpl(object):
  """ Implementation class for extracting data from VASP output """

  def __init__(self, directory = "", comm = None):
    from .. import files
    """ Initializes the extraction class. 

        @param comm: MPI group communicator. Extraction will be performed
                        for all procs in the group. In serial mode, comm can
                        be None.
        @param comm: boost.mpi.Communicator
        @param directory: path to the directory where the VASP output is located.
        @type directory: str
    """
    self.directory = directory
    """ Directory where to check for output. """
    self.comm = comm
    """ MPI group communicator. """
    self.OUTCAR = files.OUTCAR
    """ Filename of the OUTCAR file from VASP. """
    self.CONTCAR = files.CONTCAR
    """ Filename of the CONTCAR file from VASP. """
    self.FUNCCAR = files.FUNCCAR
    """ Filename of the FUNCCAR file containing the pickled functional. """
    
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

  def uncache(self): 
    """ Removes cached results.

        After this outputs are re-read from file.
    """
    from ...opt.decorators import uncache
    uncache(self)


  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_energy_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    from os.path import exists, join
    from re import compile, X as re_X
    import quantities as pq

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                           energy\(sigma->0\)\s+=\s+(\S+)""", re_X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(2))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result * pq.eV

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_energy(self):
    """ Greps total energy from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}."""
    from os.path import exists, join
    from re import compile, X as re_X
    import quantities as pq

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""energy\s+without\s+entropy\s*=\s*(\S+)\s+
                           energy\(sigma->0\)\s+=\s+(\S+)""", re_X)
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))
    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result * pq.eV

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_free_energy(self):
    """ Greps total free energy from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    from os.path import exists, join
    from re import compile
    import quantities as pq

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""free\s+energy\s+TOTEN\s*=\s*(\S+)\s+eV""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result * pq.eV

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_fermi_energy(self):
    """ Greps fermi energy from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    from os.path import exists, join
    from re import compile
    import quantities as pq

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      regex = compile( r"""E-fermi\s*:\s*(\S+)""" )
      for line in file:
        match = regex.search(line)
        if match != None: result = float(match.group(1))

    if result == None: raise RuntimeError, "File %s is incomplete.\n" % (path)
    return result * pq.eV


  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_name(self):
    """ Gets name of system from OUTCAR. """
    from os.path import exists, join
    from re import compile

    species_in = self.species

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)
    
    with open(path, "r") as file: 
      line_re = compile("^\s*POSCAR\s*=\s*")
      for line in file: 
        if line_re.search(line) != None: break
      return line[line.find("=")+1:].rstrip().strip()


  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_structure_data(self):
    """ Greps cell and positions from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
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

  @make_cached
  def _get_structure(self):
    """ Greps structure and total energy from OUTCAR. """
    from os.path import exists, join
    from re import compile
    from numpy import array, zeros
    from quantities import eV
    from ...crystal import Structure

    species_in = self.species
    cell, atoms = self._get_structure_data()

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

  @make_cached
  def _get_contcar_structure(self):
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

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_ions_per_type(self):
    """ Greps species from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    from os.path import exists, join
    from re import compile, X as re_X

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = []
    with open(path, "r") as file:
      regex = compile(r"""\s*ions\s+per\s+type\s*=""")
      for line in file:
        match = regex.search(line)
        if match == None: continue
        return [int(u) for u in line.split()[4:]]
    return None

  

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_species(self):
    """ Greps species from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    from os.path import exists, join
    from re import compile, X as re_X

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
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

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_fft(self):
    """ Greps recommended or actual fft setting from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
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


      
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_kpoints(self):
    """ Greps k-points from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. 
    
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

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_multiplicity(self):
    """ Greps multiplicity of each k-point from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
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


  def _get_eigocc(self,which):
    """ Implementation of _get_eigenvalues and _get_occupations """
    import re 
    from os.path import exists, join
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
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

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_eigenvalues(self):
    """ Greps eigenvalues of each band and kpoint from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}.

        Returns a two-dimension numpy nxm array of eigenvalues, with n the
        number of kpoints and m the number of bands.
    """
    import quantities as pq
    return self._get_eigocc(1) * pq.eV

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_occupations(self):
    """ Greps occupations according to k-point and\
        band index from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}.

        Returns a two-dimension numpy nxm array of occupations, with n the
        number of kpoints and m the number of bands.
    """
    import quantities as pq
    return self._get_eigocc(2) * pq.elementary_charge

  def _get_pressures(self, which):
    """ Greps pressure from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} """
    import re 
    from os.path import exists, join
    import quantities as pq

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    result = None
    with open(path, "r") as file:
      found = re.compile(r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+"""
                          """Pullay\s+stress\s*=\s*(\S+)\s*kB""", re.X )
      for line in file:
        match = found.search(line)
        if match != None: result = float(match.group(which))
    return result * pq.kbar

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_pressure(self):
    """ Greps pressure from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. """
    return self._get_pressures(1)

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_pulay_pressure(self):
    """ Greps pulay pressure from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} """
    return self._get_pressures(2)

  def _get_partial_charges_magnetization(self, grep):
    """ Greps partial charges from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} 

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
    print index, len(lines)
    if index == len(lines): return None
    index -= 4
    line_re = re.compile(r"""^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$""")
    for i in xrange(0, index): 
      match = line_re.match(lines[-index+i])
      if match == None: break
      result.append([float(match.group(j)) for j in range(1, 5)])
    return array(result, dtype="float64")


  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_partial_charges(self):
    """ Greps partial charges from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} 

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    return self._get_partial_charges_magnetization(r"""\s*total\s+charge\s*$""")

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_magnetization(self):
    """ Greps partial charges from L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} 

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    return self._get_partial_charges_magnetization(r"""^\s*magnetization\s*\(x\)\s*$""")

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_moment(self):
    """ Returns magnetic moment from OUTCAR. """
    from os.path import exists, join
    from re import compile
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    with open(path, "r") as file: lines = file.readlines()
    found = compile(r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$""") 
    for line in lines[::-1]:
      match = found.match(line)
      if match != None: return float(match.group(2))


  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_nb_electrons(self):
    """ Returns number of electrons from OUTCAR. """
    from os.path import exists, join
    from re import compile
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError, "File %s does not exist.\n" % (path)

    with open(path, "r") as file: lines = file.readlines()
    found = compile(r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$""") 
    for line in lines[::-1]:
      match = found.match(line)
      if match != None: return float(match.group(1))
    
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_success(self):
    """ Checks that VASP run has completed. 

        At this point, checks for the existence of
        L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>} and
        L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}. Then checks that timing stuff
        is present at end of L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}.
    """
    from os.path import exists, join
    import re

    for path in [self.OUTCAR]:
      if self.directory != "": path = join(self.directory, path)
      if not exists(path): return False
      
    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)

    with open(path, "r") as file:
      regex = re.compile(r"""General\s+timing\s+and\s+accounting
                             \s+informations\s+for\s+this\s+job""", re.X)
      for line in file:
        if regex.search(line) != None: return True
    return False

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_functional(self):
    """ Returns vasp functional used for calculation.

        Requires file L{FUNCCAR} to be present.
    """
    from os.path import exists, join
    from cPickle import load
    path = self.FUNCCAR if len(self.directory) == 0 else join(self.directory, self.FUNCCAR)
    if not exists(path): return None
    with open(path, "r") as file: return load(file)

  @make_cached
  @broadcast_result(attr=True, which=0)
  def _get_electropot(self):
    """ Average atomic electrostatic potentials.
    
        Requires file L{OUTCAR <lada.vasp.extract._ExtractImpl.OUTCAR>}.
        :return: an array with an entry for each atom. 
    """
    from os.path import exists, join
    import re
    from numpy import array

    path = self.OUTCAR if len(self.directory) == 0 else join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError("File %s does not exist.\n" % (path))

    with open(path, "r") as file:
      regex = re.compile(r"""average\s+\(electrostatic\)\s+potential\s+at\s+core""", re.X)
      for line in file:
        if regex.search(line) != None: break
      try:
        file.next()
        file.next()
      except:
        print "No average potential."
        return None
      result = []
      for line in file:
        data = line.split()
        if len(data) == 0: break
        result.extend( [float(u) for i, u in enumerate(data) if i % 2 == 1] )
        
    return array(result, dtype="float64")


  def __getstate__(self):
    from os.path import relpath
    d = self.__dict__.copy()
    if "comm" in d: del d["comm"]
    if "directory" in d: d["directory"] = relpath(d["directory"])
    return d
  def __setstate__(self, arg):
    self.__dict__.update(arg)
    self.comm = None

