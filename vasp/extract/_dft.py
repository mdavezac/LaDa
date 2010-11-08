""" Subpackage containing extraction methods for VASP-DFT data from output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ...opt.decorators import make_cached, broadcast_result

class Extract(object):
  """ Implementation class for extracting data from VASP output """

  def __init__(self, directory = None, comm = None):
    """ Initializes the extraction class. """
    super(Extract, self).__init__()
    
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
  def energies_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    from numpy import array
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    try: result = [float(u.group(2)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find energies in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find energy in OUTCAR")
    return array(result) * eV

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

  @property
  def energy(self): 
    """ Alias for total_energy. """
    return self.total_energy

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
  def pressures(self):
    """ Greps all pressures from OUTCAR """
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    try: result = [float(u.group(1)) for u in self._search_OUTCAR(regex)]
    except TypeError: raise RuntimeError("Could not find pressures in OUTCAR")
    assert len(result) != 0, RuntimeError("Could not find pressures in OUTCAR")
    return result * kB

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

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def alphabet(self):
    """ Greps alpha+bet from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find alpha+bet in OUTCAR")
    return float(result.group(3))

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def xc_g0(self):
    """ Greps XC(G=0) from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result != None, RuntimeError("Could not find xc(G=0) in OUTCAR")
    return float(result.group(2))

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
