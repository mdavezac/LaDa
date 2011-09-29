""" Extracts VASP-GW output """
__docformat__  = 'restructuredtext en'
__all__ = ['Extract']
from ...opt.decorators import make_cached
from ...opt.json import array as json_array,\
                        array_with_unit as json_array_with_unit,\
                        section as json_section
from quantities import eV

class Extract(object):
  """ Implementation of GW extractor. """

  def __init__(self): 
    """ Initializes the extraction class. """
    object.__init__(self)

  @property
  @json_array_with_unit("float64", eV)
  def eigenvalues(self):
    """ Deprecated. """
    from warnings import warn
    warn( DeprecationWarning('eigenvalues attribute is deprecated in favor of qp_eigenvalues.'), \
          stacklevel=2 )
    return self.qp_eigenvalues

  @property
  @json_section("output")
  @json_array_with_unit("float64", eV)
  @make_cached
  def dft_eigenvalues(self):
    """ Greps DFT eigenvalues from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(1), dtype="float64") * eV

  @property
  @json_section("output")
  @json_array_with_unit("float64", eV)
  @make_cached
  def qp_eigenvalues(self):
    """ Greps quasi-particle eigenvalues from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(2), dtype="float64") * eV

  @property
  @json_section("output")
  @json_array_with_unit("float64", eV)
  @make_cached
  def self_energies(self):
    """ Greps self-energies of each eigenvalue from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64") * eV
    return array(self._unpolarized_values(3), dtype="float64") * eV

  @property
  @json_section("output")
  @json_array("float64")
  @make_cached
  def occupations(self):
    """ Greps occupations from OUTCAR.

        In spin-polarized cases, the leading dimension of the numpy array are
        spins, followed by kpoints, and finally with bands. In spin-unpolarized
        cases, the leading dimension are the kpoints, followed by the bands.
    """
    from numpy import array
    if self.ispin == 2: return array(self._spin_polarized_values(1), dtype="float64")
    return array(self._unpolarized_values(3), dtype="float64")
