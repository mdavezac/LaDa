"""" Creates diagram of formation enthalpies with respect to Fermi energy. """
from lada.jobs import AbstractMassExtract

class Enthalpy(object):
  """ Enthalpy of a single defect. """
  def __init__(self, extract, epsilon = 1e0, host = None, pa_maxdiff=-8):
    """ Initializes an enthalpy function. """
    # extraction object.
    self.extract = extract.copy(unix_re=False)
    if self.extract.excludes == None: self.extract.excludes = [".*relax_*"]
    else: self.extract.excludes.append(".*relax_*$")

    self.epsilon = epsilon
    """ Dimensionless dielectric constant. """
    self._host = host
    """ Host extraction object. 

        If None, will determine it from the directory structure.
    """
    self.pa_maxdiff = pa_maxdiff
    """ Potential alignment parameter. """

  @property
  def rootdir(self):
    """ Root directory of defects. """
    return self.extract.rootdir

  def _charge_correction(self, extract):
    """ Returns the charge correction. """
    assert len(extract.jobs) == 1, ValueError("extract should return naked objects.")
    assert extract.naked_end,      ValueError("extract should return naked objects.")
    return extract.charge_corrections / self.epsilon 

  def _potential_alignment(self, extract):
    """ Returns the charge correction. """
    from lada.crystal.point_defects import potential_alignment
    assert len(extract.jobs) == 1, ValueError("extract should return naked objects.")
    assert extract.naked_end,      ValueError("extract should return naked objects.")
    return potential_alignment(extract, self.host, self.pa_maxdiff)

  def _band_filling(self, extract):
    from lada.crystal.point_defects import band_filling
    assert len(extract.jobs) == 1, ValueError("extract should return naked objects.")
    assert extract.naked_end,      ValueError("extract should return naked objects.")
    return band_filling(extract, self.host.cbm + self._potential_alignment(extract))

  def _corrected(self, extract):
    """ Corrected formation enthalpy. """
    print extract.view
    print extract.total_energy, self.host.total_energy, self._charge_correction(extract),\
          self._potential_alignment(extract), self._band_filling(extract)
    return   extract.total_energy\
           - self.host.total_energy \
           + self._charge_correction(extract)\
           + self._potential_alignment(extract)\
           + self._band_filling(extract)


  def _charged_states(self):
    """ Yields extraction routine toward each charge states. 

        If a charge state has children, then only the lowest energy calculation
        is returned.
    """
    from operator import itemgetter
    for child in self.extract.children: # Each child is a different charge state.
      child.naked_end = False
      lowest = sorted(child.total_energies.items(), key=itemgetter(1))[0][0]
      result = child[lowest]
      result.naked_end = True
      yield result

  @property
  def host(self):
    """ Returns extraction object towards the host. """
    if self._host == None: 
      host = self.extract['../..' if self._is_site == None else '../../..']
      host = self.copy(excludes=[".*PointDefects"], naked_end=False)
      host.excludes.extend(host.excludes)
      lowest = sorted(child.total_energies.items(), key=itemgetter(1))[0][0]
      self._host = [u for u in self.extract[lowest].values()]
      assert len(self._host) == 1
      self._host = self._host[0]
    return self._host

  @property 
  def _site(self):
    """ Returns site number or None. """
    from re import match
    regex = match("^site_(\d+)$", self.extract.view.split()[-1])
    return int(regex.group(1)) if regex != None else None

  @property 
  def name(self):
    """ Name of the defect. """
    return self.extract.view.split('/')[-2 if self._site != None else -1]

  @property
  def is_vacancy(self):
    """ True if this is a vacancy. """
    from re import match
    return match("vacancy_[A-Z][a-z]?", self.name) != None

  @property
  def is_interstitial(self):
    """ True if this is an interstitial. """
    from re import match
    return match("[A-Z][a-z]?_interstitial_\S+", self.name) != None

  @property
  def is_substitution(self):
    """ True if this is a substitution. """
    from re import match
    return match("[A-Z][a-z]?_on_[A-Z][a-z]?", self.name) != None

  @property
  def species(self):
    """ List of species involved in this defect. """
    from re import match
    if self.is_vacancy:
      return [match("vacancy_([A-Z][a-z])?", self.name).group(1)]
    elif self.is_interstitial:
      return [match("[A-Z][a-z]?_interstitial_(\S+)", self.name).group(1)]
    else: 
      found = match("([A-Z][a-z])?_on_([A-Z][a-z])?", self.name)
      return [match.group(1), match.group(2)]

  def _lines(self):
    """ Returns lines composed by the different charge states. """
    from numpy import array
    from quantities import elementary_charge as e
    lines = []
    states = set()
    for state in self._charged_states():
      assert state.charge not in states,\
             RuntimeError("Found more than one calculation for the same charge state.")
      states.add(state.charge)
      lines.append(array([self._corrected(state), state.charge]))
    return lines

  def _all_intersections(self, _lines):
    """ Returns all intersection points between vbm and cbm, ordered. """
    from numpy import array
    from quantities import eV
    vbm = self.host.vbm
    cbm = self.host.cbm
    result = []
    for i, (b0, a0) in enumerate(_lines[:-1]):
      for b1, a1 in _lines[i:]: result.append( (b0 - b1) / (a1 - a0) )
    result = [u for u in sorted(result) if u + 1e-6 * eV > vbm]
    result = [u for u in result if u - 1e-6 * eV > cbm]
    result.append(cbm)
    result.insert(0, vbm)
    return array([u for u in result]).rescaled(eV) 

  def lines(self):
    """ Lines forming the formation enthalpy diagram. 
    
        :return: A list of 2-tuples with the first item b and the second a (a*x+b).
    """
    from numpy import array
    from quantities import eV
    _lines = self._lines()
    intersections = self._all_intersections(_lines)
    # adds line before vbm
    func  = lambda x: x[0] + (intersections[0]-eV)*x[1] 
    lines = [ min(_lines, key=func)  ]

    # now look for lines up to cbm
    for i, intersection in enumerate(intersections[1:]):
      func  = lambda x: x[0] + (intersection-intersections[i])*x[1] 
      min_line = min(lines, key=func)
      if    abs(min_line[0] - lines[-1][0]) > 1e-12*eV \
         or abs(min_line[1] - lines[-1][1]) > 1e-12:
        lines.append([min_line[0].rescale(eV). min_line[1]])

    # adds line after cbm
    func  = lambda x: x[0] + (intersections[-1]+eV)*x[1] 
    min_line = min(_lines, key=func)
    if    abs(min_line[0] - lines[-1][0]) > 1e-12*eV \
       or abs(min_line[1] - lines[-1][1]) > 1e-12:
      lines.append([min_line[0].rescale(eV). min_line[1]])
    return lines

  def __call__(self, fermi, mu = None):
    """ Returns formation enthalpy for given fermi energy and mu. """
    from quantities import eV
    if mu == None: mu = 0e0
    elif self.is_substitutional: mu = mu[self.species[0]] - mu[self.specie[1]]
    else: mu = mu[self.species[0]]
    if hasattr(mu, "rescale"): mu = mu.rescale(eV)
    else: mu = mu * eV
    if hasattr(fermi, "rescale"): fermi = fermi.rescale(eV)
    else: fermi = fermi * eV
    return (min(x[0]+fermi*x[1] for x in self.lines()) + self.n * mu).rescale(eV)

  @property
  def latex_label(self):
    """ A label in LaTex format. """
    from re import match
    if self.is_interstitial:
      found = match("([A-Z][a-z]?)_interstitial_(.+)$", self.name) 
      return r"{0}$_{{ \\mathrm{{ {1} }} }}$".format(found.group(1), found.group(2))
    if self.is_substitutional:
      found = match("([A-Z][a-z]?)_on_([A-Z][a-z]?)", self.name) 
      site = self._is_site
      if site == None:
        return r"{0}$_{{ \\mathrm{{ {1} }} }}$".format(found.group(1), found.group(2))
      else:
        return r"{0}$_{{ \\mathrm{{ {1} }}_{{ {2} }} }}$"\
               .format(found.group(1), found.group(2), site)


class Enthalpies(AbstractMassExtract):
  """ Enthalpy for a series of defects for a given material and lattice. """
  def __init__(self, path=None, epsilon = 1e0, pa_maxdiff=0.5, **kwargs):
    """ Initializes an enthalpy function. """
    super(Enthalpies, self).__init__(**kwargs)

    from lada.vasp import MassExtract
    self.epsilon = epsilon
    """ Dimensionless dielectric constant. """
    self.pa_maxdiff = pa_maxdiff
    """ Potential alignment parameter. """

    self.massextract = MassExtract(path, unix_re=False, excludes=[".*relax_*"])
    """ Mass extraction object from which all results are pulled. """
    
    self.host = self._get_host()
    """ Result of host calculations. """

  @property
  def rootdir(self): 
    """ Path to the root-directory containing the poin-defects. """
    return self.massextract.rootdir
  @rootdir.setter
  def rootdir(self, value): self.massextract.rootdir = value
    
  def _get_host(self):
    """ Returns extraction object towards the host. """
    from operator import itemgetter
    host = self.massextract.copy(excludes=[".*PointDefects"])
    host.excludes.extend(self.massextract.excludes)
    lowest = sorted(host.total_energies.items(), key=itemgetter(1))[0][0]
    host = [u for u in host[lowest].values()]
    assert len(host) == 1
    return host[0]

  def walk_through(self):
    """ Walks through point-defects only. """
    for child in self.massextract["PointDefects"].children:
      # looks for site_n
      if len(child["site_\d+"].jobs) != 0:
        assert len(child["site_\d+"].jobs) == len(child.jobs),\
               RuntimeError("Don't understand directory structure of {0}.".format(child.view))
        for site in child.children: # should site specific defects.
          result = Enthalpy(child, self.epsilon, self.host, self.pa_maxdiff)
          # checks this is a point-defect.
          if result.is_interstitial or result.is_vacancy or result.is_substitution:
            yield site.view, result
      else:
        result = Enthalpy(child, host=self.host, pa_maxdiff=self.pa_maxdiff, epsilon = self.epsilon)
        # checks if this is a valid point-defect.
        if result.is_interstitial or result.is_vacancy or result.is_substitution:
          yield child.view, result

  def __call__(self, fermi, mu=None):
    """ Returns dictionary of point-defect formation enthalpies. """
    results = {}
    for name, defect in self._regex_extractors(): results[name] = defect(fermi, mu)
    return results
      

  def plot_enthalpies(self, mu=None, **kwargs):
    """ Plots diagrams using matplotlib. """
    from quantities import eV
    try: import matplotlib.pyplot as plt
    except ImportError: 
      print "No matplotlib module."
      return
    from operator import itemgetter

    # finds limits of figure
    xlim = float(self.host.vbm.rescale(eV)), float(self.host.cbm.rescale(eV)) 
    ylim = min(self(xlim[0], mu).items(), key=itemgetter(1)),\
           max(self(xlim[0], mu).items(), key=itemgetter(1) )
    ylim = float(ylim[0].rescale(eV)), float(ylim[1].rescale(eV))
    # creates figures and axes.
    figure = plt.figure()
    figure.add_subplot(111, xlim=(self.vbm, self.cbm), ylim=ylim)

    # loop over defects.
    for name, defect in self:
      # finds intersection points. 
      x = [self.host.vbm.rescale(eV) - 5e0*eV]
      lines = defect.lines()
      for i in range(len(lines)-1):
        (b0, a0), (b1, a1) = lines[i], lines[i+1]
        intersections.append( ((b0-b1)/(a1-a0)).rescale(eV) )
      x = [self.host.cbm.rescale(eV) + 5e0*eV]

      # Now draws lines. 
      lines.append(lines[-1])
      y = [u[0] + u[1] * xx for u, xx in zip(lines, x)]
      figure.plot(x, y, label=defect.latex_label, **kwargs)
      
        









