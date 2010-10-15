"""" Creates diagram of formation enthalpies with respect to Fermi energy. """

from lada.vasp import MassExtract
class Enthalpy(MassExtract):
  """ Enthalpy of a single defect. """
  def __init__(self, path, epsilon = 1e0, host = None, pa_maxdiff=0.5, **kwargs):
    """ Initializes an enthalpy function. """
    unix_re = kwargs.pop("unix_re", False)
    assert not unix_re, ValueError("Enthalpy only works with python regex.")
    super(Enthalpy, self).__init__(path, unix_re=False, **kwargs)

    # excludes relaxation stuff.
    if self.excludes == None: self.excludes = [".*relax_*"]
    else: self.excludes.append(".*relax_*$")

    self.epsilon = epsilon
    """ Dimensionless dielectric constant. """
    self._host = host
    """ Host extraction object. 

        If None, will determine it from the directory structure.
    """
    self.pa_maxdiff = pa_maxdiff
    """ Potential alignment parameter. """

  @property 
  def vbm(self):
    """ vbm of the host. """
    from numpy import max
    return max(self.host.eigenvalues[:, valence])

  @property
  def cbm(self):
    """ cbm of the host. """
    from numpy import min
    return min(self.host.eigenvalues[:, valence+1])

  @property
  def _charge_correction(self):
    """ Returns the charge correction. """
    from lada.crystal.point_defect import charge_correction
    assert len(self.jobs) == 0, RuntimeError("Cannot use this property.")
    return charge_correction(self.cell, charge=self.charge, self.epsilon)

  @property
  def _potential_alignment(self):
    """ Returns the charge correction. """
    from lada.crystal.point_defect import potential_alignment
    assert len(self.jobs) == 0, RuntimeError("Cannot use this property.")
    return potential_alignment(self, self.host, self.pa_maxdiff)

  @property
  def _band_filling(self):
    from lada.crystal.point_defect import band_filling
    from pq import elementary_charge as e
    from numpy import max
    assert len(self.jobs) == 0, RuntimeError("Cannot use this property.")
    valence = int(self.valence.rescale(e)+1e-12)
    return band_filling(extract, self.cbm + self._potential_aligment(extract))

  @property
  def _corrected(self):
    """ Corrected formation enthalpy. """
    assert len(self.jobs) == 0, RuntimeError("Cannot use this property.")
    return self.total_energy - self.host.total_energy \
           + self._charge_correction + self._potential_alignemnt + self._band_filling


  def _charged_states(self):
    """ Yields extraction routine toward each charge states. 

        If a charge state has children, then only the lowest energy calculation
        is returned.
    """
    for child in self.children: # Each child is a different charge state.
      child.naked_end = False
      lowest = sorted(child.total_energies.items(), key=itemgetter(1))[0][0]
      result = child[lowest]
      result.naked_end = True
      yield result

  @property
  def host(self):
    """ Returns extraction object towards the host. """
    if self._host == None: 
      from os.path import dirname
      from copy import deepcopy
      path = dirname(self.path)
      if self._is_site != None: path = dirname(self.path)
      assert path.split('/')[-1] == "PointDefects", RuntimeError("Could find host material.")
      path = dirname(path)
      host = self.copy(path=path)
      host.excludes = deepcopy(host.excludes)
      host.excludes.append(".*PointDefects")
      lowest = sorted(child.total_energies.items(), key=itemgetter(1))[0][0]
      self._host = host[lowest]
      self._host.naked_end = True
    return self._host

  @property 
  def _site(self):
    """ Returns site number or None. """
    from re import match
    regex = match("^site_(\d+)$", self.path.split()[-1])
    return int(regex.group(1)) if regex != None else None

  @property 
  def name(self):
    """ Name of the defect. """
    from os.path import relative, dirname
    result = relative(self.path, dirname(self.path))
    return result if self._site == None else relative(result, dirname(self.path))

  @property
  def is_vacancy(self):
    """ True if this is a vacancy. """
    from re import match
    return match("Vacancy_[A-Z][a-z]?", self.name) != None

  @property
  def is_insterstitial(self):
    """ True if this is an interstitial. """
    from re import match
    return match("[A-Z][a-z]?_interstitial_\S+", self.name) != None

  @property
  def is_substition(self):
    """ True if this is a substitution. """
    from re import match
    return match("[A-Z][a-z]?_on_[A-Z][a-z]?", self.name) != None

  @property
  def species(self):
    """ List of species involved in this defect. """
    from re import match
    if self.is_vacancy:
      return [match("Vacancy_([A-Z][a-z])?", self.name).group(1)]
    elif self.is_insterstitial:
      return [match("[A-Z][a-z]?_interstitial_(\S+)", self.name).group(1)]
    else: 
      found = match("([A-Z][a-z])?_on_([A-Z][a-z])?", self.name)
      return [match.group(1), match.group(2)]

  def _lines():
    """ Returns lines composed by the different charge states. """
    from quantitie import elementary_charge as e
    lines = []
    states = set()
    for state in self._charged_state():
      assert state.charge not in states,\
             RuntimeError("Found more than one calculation for the same charge state.")
      states.add(state.charge)
      lines.append((state._corrected, float(state.charge.rescale(e))))

  def _all_intersections(self):
    """ Returns all intersection point, ordered. """
    from quantities import eV
    vbm = self.vbm
    cbm = self.cbm
    intersections = [cbm, vbm]
    lines = self._lines()
    for i, (b0, a0) in enumerate(lines[:-1]):
      for b1, a1 in lines[i:]: intersections.append( (b0 - b1) / (a1 - a0) )
    intersections =  sorted(intersections)
    for i, s in enumerate(intersections):
      if float(abs(s-vbm).simplified) < 1e-12: break
    for j, s in enumerate(intersections[i:]):
      if float(abs(s-cbm).simplified) < 1e-12: break
    return array([u for u in intersections[i:j+1].rescaled(eV))]) 

  def lines(self):
    """ Lines forming the formation enthalpy diagram. 
    
        :return: A list of 2-tuples with the first item b and the second a (a*x+b).
    """
    lines = None
    for intersection in intersections:
      min_line = min(lines, key=lambda x: x[0] + intersection*x[1])
      if lines == None: lines = [min_line]
      elif any(abs(min_line - lines) > 1e-12): lines.append(min_line)
    return lines

  def DeltaH(fermi, mu = None):
    """ Returns formation enthalpy for given fermi energy and mu. """
    from quantities import eV
    if mu == None: mu = 0e0
    elif self.is_substitutional: mu = mu[self.species[0]] - mu[self.specie[1]]
    else: mu = mu[self.species[0]]
    return min([ x[0]+fermi*x[1] for x in self.lines ]) + self.n * mu * eV


 


