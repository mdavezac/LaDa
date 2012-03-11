""" Subpackage containing extraction methods for VASP-DFT data from output. """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ._common import Extract as ExtractCommon
from ...functools import make_cached

class Extract(ExtractCommon):
  """ Implementation class for extracting data from VASP output """

  def __init__(self):
    """ Initializes the extraction class. """
    super(Extract, self).__init__()
    
  @property
  @make_cached
  def energy_sigma0(self):
    """ Greps total energy extrapolated to $\sigma=0$ from OUTCAR. """
    from quantities import eV
    regex = """energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find sigma0 energy in OUTCAR")
    return float(result.group(2)) * eV

  @property
  @make_cached
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
  def fermi0K(self):
    """ Fermi energy at zero kelvin.
    
        This procedure recomputes the fermi energy using a step-function. 
        It avoids negative occupation numbers. As such, it may be different
        from the fermi energy given by vasp, depeneding on the smearing and the
        smearing function.
    """
    from operator import itemgetter
    from numpy import rollaxis, sum

    if self.ispin == 2: eigs = rollaxis(self.eigenvalues, 0, 2)
    else: eigs = self.eigenvalues
    array = [(e, m) for u, m in zip(eigs, self.multiplicity) for e in u.flat]
    array = sorted(array, key=itemgetter(0))
    occ = 0e0
    factor = (2.0 if self.ispin == 1 else 1.0)/float(sum(self.multiplicity))
    for i, (e, w) in enumerate(array):
      occ += w*factor
      if occ >= self.valence - 1e-8: return e * self.eigenvalues.units
    return None

  @property
  def halfmetallic(self):
    """ True if the material is half-metallic. """
    from numpy import max, abs
    if self.ispin == 1: return False

    if self.cbm - self.vbm > 0.01 * self.cbm.units: return False

    fermi = self.fermi0K
    vbm0 = max(self.eigenvalues[0][self.eigenvalues[0]<=float(fermi)+1e-8]) 
    vbm1 = max(self.eigenvalues[1][self.eigenvalues[1]<=float(fermi)+1e-8])
    return abs(float(vbm0-vbm1)) > 2e0 * min(float(fermi - vbm0), float(fermi - vbm1))

  @property
  def cbm(self):
    """ Returns Condunction Band Minimum. """
    from numpy import min
    return min(self.eigenvalues[self.eigenvalues>self.fermi0K+1e-8*self.fermi0K.units])

  @property
  def vbm(self):
    """ Returns Valence Band Maximum. """
    from numpy import max
    from quantities import eV
    return max(self.eigenvalues[self.eigenvalues<=self.fermi0K+1e-8*self.fermi0K.units])

  @property
  @make_cached
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
  def total_energy(self):
    """ Greps total energy from OUTCAR."""
    from quantities import eV
    regex = """energy\s+without\s+entropy=\s*(\S+)\s+energy\(sigma->0\)\s+=\s+(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find energy in OUTCAR")
    return float(result.group(1)) * eV

  energy = total_energy
  """ Alias for total_energy. """

  @property
  @make_cached
  def fermi_energy(self):
    """ Greps fermi energy from OUTCAR. """
    from quantities import eV
    regex = r"""E-fermi\s*:\s*(\S+)"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find fermi energy in OUTCAR")
    return float(result.group(1)) * eV

  @property
  @make_cached
  def moment(self):
    """ Returns magnetic moment from OUTCAR. """
    regex = r"""^\s*number\s+of\s+electron\s+(\S+)\s+magnetization\s+(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find magnetic moment in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
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
  def pressure(self):
    """ Greps last pressure from OUTCAR """
    from quantities import kbar as kB
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find pressure in OUTCAR")
    return float(result.group(1)) * kB

  @property
  @make_cached
  def alphabet(self):
    """ Greps alpha+bet from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find alpha+bet in OUTCAR")
    return float(result.group(3))

  @property
  @make_cached
  def xc_g0(self):
    """ Greps XC(G=0) from OUTCAR """
    regex = r"""^\s*E-fermi\s*:\s*(\S+)\s+XC\(G=0\)\s*:\s*(\S+)\s+alpha\+bet\s*:(\S+)\s*$"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find xc(G=0) in OUTCAR")
    return float(result.group(2))

  @property
  @make_cached
  def pulay_pressure(self):
    from quantities import kbar as kB
    """ Greps pressure from OUTCAR """
    regex = r"""external\s+pressure\s*=\s*(\S+)\s*kB\s+Pullay\s+stress\s*=\s*(\S+)\s*kB"""
    result = self._find_last_OUTCAR(regex) 
    assert result is not None, RuntimeError("Could not find pulay pressure in OUTCAR")
    return float(result.group(2)) * kB

  @property
  @make_cached
  def fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    from re import search, M as re_M
    from numpy import array
    regex = r"""(?!I would recommend the setting:\s*\n)"""\
             """\s*dimension x,y,z NGX =\s+(\d+)\s+NGY =\s+(\d+)\s+NGZ =\s+(\d+)"""
    with self.__outcar__() as file: result = search(regex, file.read(), re_M)
    assert result is not None, RuntimeError("Could not FFT grid in OUTCAR.""")
    return array([int(result.group(1)), int(result.group(2)), int(result.group(3))])

  @property
  @make_cached
  def recommended_fft(self):
    """ Greps recommended or actual fft setting from OUTCAR. """
    from re import search, M as re_M
    from numpy import array
    regex = r"""I would recommend the setting:\s*\n"""\
             """\s*dimension x,y,z NGX =\s+(\d+)\s+NGY =\s+(\d+)\s+NGZ =\s+(\d+)"""
    with self.__outcar__() as file: result = search(regex, file.read(), re_M)
    assert result is not None, RuntimeError("Could not recommended FFT grid in OUTCAR.""")
    return array([int(result.group(1)), int(result.group(2)), int(result.group(3))])

  def _get_partial_charges_magnetization(self, grep):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included. Implementation also used for magnetization.
    """
    import re 
    from numpy import array

    result = []
    with self.__outcar__() as file: lines = file.readlines()
    found = re.compile(grep) 
    for index in xrange(1, len(lines)+1):
      if found.search(lines[-index]) is not None: break 
    if index == len(lines): return None
    index -= 4
    line_re = re.compile(r"""^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$""")
    for i in xrange(0, index): 
      match = line_re.match(lines[-index+i])
      if match is None: break
      result.append([float(match.group(j)) for j in range(1, 5)])
    return array(result, dtype="float64")

  @property
  @make_cached
  def partial_charges(self):
    """ Greps partial charges from OUTCAR.

        This is a numpy array where the first dimension is the ion (eg one row
        per ion), and the second the partial charges for each angular momentum.
        The total is not included.
    """
    return self._get_partial_charges_magnetization(r"""\s*total\s+charge\s*$""")

  @property
  @make_cached
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
  def electropot(self):
    """ Greps average atomic electrostatic potentials from OUTCAR. """
    from re import compile, X as reX
    from numpy import array
    from quantities import eV

    with self.__outcar__() as file: lines = file.readlines()
    regex = compile(r"""average\s+\(electrostatic\)\s+potential\s+at\s+core""", reX)
    for i, line in enumerate(lines[::-1]):
      if regex.search(line) is not None: break
    assert -i + 2 < len(lines), RuntimeError("Could not find average atomic potential in file.")
    regex = compile(r"""(?:\s|\d){8}\s*(-?\d+\.\d+)""")
    result = []
    for line in lines[-i+2:]:
      data = line.split()
      if len(data) == 0: break
      result.extend([m.group(1) for m in regex.finditer(line)])
        
    return array(result, dtype="float64") * eV

  @property 
  @make_cached
  def electronic_dielectric_constant(self):
    """ Electronic contribution to the dielectric constant. """
    from re import M as multline
    from numpy import array
    regex = r"\s*MACROSCOPIC\s+STATIC\s+DIELECTRIC\s+TENSOR\s*\(including local field effects in DFT\)\s*\n"\
            r"\s*-+\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*-+\s*\n"
    result = self._find_last_OUTCAR(regex, multline)
    assert result is not None, RuntimeError('Could not find dielectric tensor in output.')
    return array([result.group(i) for i in range(1,10)], dtype='float64').reshape((3,3))

  @property 
  @make_cached
  def ionic_dielectric_constant(self):
    """ Ionic contribution to the dielectric constant. """
    from re import M as multline
    from numpy import array
    regex = r"\s*MACROSCOPIC\s+STATIC\s+DIELECTRIC\s+TENSOR\s+IONIC\s+CONTRIBUTION\s*\n"\
            r"\s*-+\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*(\S+)\s+(\S+)\s+(\S+)\s*\n"\
            r"\s*-+\s*\n"
    result = self._find_last_OUTCAR(regex, multline)
    assert result is not None, RuntimeError('Could not find dielectric tensor in output.')
    return array([result.group(i) for i in range(1,10)], dtype='float64').reshape((3,3))

  @property 
  def dielectric_constant(self):
    """ Dielectric constant of the material. """
    return  self.electronic_dielectric_constant + self.ionic_dielectric_constant

  @property
  @make_cached
  def stresses(self):
    """ Returns total stress at each relaxation step. """
    from numpy import zeros, abs, dot, all, array
    from numpy.linalg import det
    from quantities import eV, J, kbar
    from re import finditer, M 
    if self.isif < 1: return None
    pattern = """\s*Total\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"""\
              """\s*in kB\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\n"""
    result = []
    with self.__outcar__() as file:
      for regex in finditer(pattern, file.read(), M): 
        stress = zeros((3,3), dtype="float64"), zeros((3,3), dtype="float64") 
        for i in xrange(2): 
          for j in xrange(3): stress[i][j, j] += float(regex.group(i*6+1+j))
          stress[i][0,1] += float(regex.group(i*6+4))
          stress[i][1,0] += float(regex.group(i*6+4))
          stress[i][1,2] += float(regex.group(i*6+4))
          stress[i][2,1] += float(regex.group(i*6+4))
          stress[i][0,2] += float(regex.group(i*6+4))
          stress[i][2,0] += float(regex.group(i*6+4))
        if sum(abs(stress[0].flatten())) > sum(abs(stress[1].flatten())):
          result.append( stress[0] * float(eV.rescale(J) * 1e22\
                         / abs(det(self.structure.cell*self.structure.scale))) * kbar)
        else: result.append(stress[1])
    return array(result)*kbar
  @property
  def stress(self):
     """ Returns final total stress. """
     return self.stresses[-1]

  @property
  @make_cached
  def forces(self):
    """ Forces on each atom. """
    from numpy import array
    from quantities import angstrom, eV
    from re import compile, finditer, M
    result = []
    pattern = """ *POSITION\s*TOTAL-FORCE\s*\(eV\/Angst\)\s*\n"""\
              """\s*-+\s*\n"""\
              """(?:\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s*\n)+"""\
              """\s*-+\s*\n"""\
              """\s*total drift"""
    with self.__outcar__() as file:
      for regex in finditer(pattern, file.read(), M): pass
    return array([u.split()[3:] for u in regex.group(0).split('\n')[2:-2]], dtype="float64")\
           * eV / angstrom

  @property
  @make_cached
  def structure(self):
    """ Greps structure from outcar. 

        Adds attributes from initial structure where possible.
        Adds magnetization and forces to atoms.
    """
    # adds attributes from initial structure.
    result = super(Extract, self).structure
    try: initial = self.initial_structure
    except: pass
    else:
      for key, value in initial.__dict__.iteritems():
        if not hasattr(result, key): setattr(result, key, value)
      for a, b in zip(result, initial):
        for key, value in b.__dict__.iteritems():
          if not hasattr(a, key): setattr(a, key, value)
    # adds magnetization.
    try: magnetization = self.magnetization
    except: pass
    else:
      if magnetization is not None:
        for atom, mag in zip(result, magnetization): atom.magmom = sum(mag)
    # adds stress.
    try: result.stress = self.stress
    except: pass
    # adds forces.
    try: forces = self.forces
    except: pass
    else:
      for force, atom in zip(forces, result): atom.force = force
    return result
