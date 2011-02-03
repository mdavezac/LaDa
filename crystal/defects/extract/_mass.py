"""" Extraction object for many point-defects, but single material. """
__docformat__ = "restructuredtext en"
__all__ = ['Material']
from ....jobs import AbstractMassExtract
from ._single import Single

class _MaterialNavigator(AbstractMassExtract):
  """ Navigates around multiple defects of a single material. """
  DefectExtractor = Single
  """ Class for extracting data from a single defect. """
  def __init__(self, path=None, epsilon = 1e0, pa_maxdiff=0.5, **kwargs):
    """ Initializes an enthalpy function. """
    from ....vasp import MassExtract as VaspMassExtract

    # possible customization of mass defect extration object.
    MassExtractor = kwargs.pop("MassExtractor", VaspMassExtract)
    """ Object type for mass extraction. """
    # possible customization of single defect extration objects.
    self.__dict__["DefectExtractor"] \
        = kwargs.pop("DefectExtractor", _MaterialNavigator.DefectExtractor)

    super(_MaterialNavigator, self).__init__(**kwargs)

    self.massextract = MassExtractor(path, unix_re=False, excludes=[".*relax_*"])
    """ Mass extraction object from which all results are pulled. """
    self.host = self._get_host()
    """ Result of host calculations. """

    # must be last. Can't use porperty setter.
    self._epsilon = epsilon
    self._pa_maxdiff = pa_maxdiff


  @property 
  def epsilon(self): 
    """ Dimensionless dielectric constant. """
    return self._epsilon
  @epsilon.setter
  def epsilon(self, value):
    self._epsilon = value 
    for v in self.itervalues(): v.epsilon = self._epsilon

  @property 
  def pa_maxdiff(self): 
    """ Dimensionless dielectric constant. """
    return self._pa_maxdiff
  @pa_maxdiff.setter
  def pa_maxdiff(self, value):
    self._pa_maxdiff = value 
    for v in self.itervalues(): v.pa_maxdiff = self._pa_maxdiff


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
    lowest = sorted(host.total_energies.iteritems(), key=itemgetter(1))[0][0]
    host = [u for u in host[lowest].itervalues()]
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
          result = self.DefectExtractor(site, self.epsilon, self.host, self.pa_maxdiff)
          # checks this is a point-defect.
          if result.is_interstitial or result.is_vacancy or result.is_substitution:
            yield site.view, result
      else:
        result = self.DefectExtractor(child, host=self.host, pa_maxdiff=self.pa_maxdiff,\
                                      epsilon = self.epsilon)
        # checks if this is a valid point-defect.
        if result.is_interstitial or result.is_vacancy or result.is_substitution:
          yield child.view, result

  def ordered_items(self):
    """ Returns items ordered by substitution, vacancy, and interstitial. """
    from operator import itemgetter
    interstitials = (u for u in self.iteritems() if u[1].is_interstitial)
    substitution  = (u for u in self.iteritems() if u[1].is_substitution)
    vacancy       = (u for u in self.iteritems() if u[1].is_vacancy)
    result = sorted(substitution, key = itemgetter(0)) 
    result.extend(sorted(vacancy, key = itemgetter(0)))
    result.extend(sorted(interstitials, key = itemgetter(0)))
    return result
  def ordered_keys(self):
    """ Returns keys ordered by substitution, vacancy, and interstitial. """
    return [u[0] for u in self.ordered_items()]
  def ordered_values(self):
    """ Returns values ordered by substitution, vacancy, and interstitial. """
    return [u[1] for u in self.ordered_items()]

        


class Material(_MaterialNavigator):
  """ Extracts data for a whole material. """
  def __init__(self, *args, **kwargs):
    """ Initializes an enthalpy function. """
    super(Material, self).__init__(*args, **kwargs)

  @property
  def cbm(self):
    """ Conduction band minimum of the host. """
    return self.host.cbm
  @property
  def vbm(self):
    """ Valence band maximum of the host. """
    return self.host.vbm

  def enthalpies(self, fermi, mu=None):
    """ Dictionary of point-defect formation enthalpies. 
    
        :Parameters:
          fermi  
            Fermi energy with respect to the host's VBM. If without
            units, assumes eV. 
          mu : dictionary or None
            Dictionary of chemical potentials. If without units, assumes eV.
            If None, chemical potential part of the formation enthalpy is
            assumed zero.

        :return: Dictionary where keys are the name of the defects, and the
          values the formation enthalpy.
    """
    from quantities import eV
    results = {}
    for name, defect in self.iteritems():
      results[name] = defect.enthalpy(fermi, mu).rescale(eV)
    return results
  
  def __str__(self): 
    """ Prints out all energies and corrections. """
    return "".join( str(value) for value in self.ordered_values() )
      
