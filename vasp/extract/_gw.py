""" Extracts VASP-GW output """
__docformat__  = 'restructuredtext en'
__all__ = ['ExtratGW']
from ._dft import _ExtractImpl
from ...opt.decorators import make_cached, broadcast_result
class ExtractGW(_ExtractImpl):
  """ Implementation of GW extractor. """

  def __init__(self, directory = None, comm = None, **kwargs):
    """ Initializes the extraction class. 

        :Parameters:
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory.
          comm : boost.mpi.communicator or None
            MPI group communicator. Extraction will be performed for all procs
            in the group. In serial mode, comm can be None.
    """
    super(ExtractGW, self).__init__(directory=directory, comm=comm, **kwargs)

  @property
  def is_dft(self): return False

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
