""" Point Charge Model functional. """
__docformat__ = "restructuredtext en"

from _pcm import Charges, LJBond, LJBonds, Clj, bond_name, minimize

__all__ = [ "Charges", "LJBond", "LJBonds", "Clj", "bond_name", "minimize" ]

def _get_gcutoff(self):
  """ Cutoff in G-space for the Ewald sum. 

      By default, units are in eV, but will adapt to whatever units from the
      package `quantities <http://packages.python.org/quantities/user/tutorial.html>`_.
  """
  from quantities import eV
  if not hasattr(self, "_cutoff_units"):
    self._cutoff_units = eV
  return (self._cutoff * eV).rescale(self._cutoff_units)
def _set_gcutoff(self, value):
  from quantities import eV
  if hasattr(value, "rescale"): 
    self._cutoff = float(value.rescale(eV))
    self._cutoff_units = value.units
  elif not hasattr(self, "_cutoff_units"):
    self._cutoff = float(value)
    self._cutoff_units = eV
  else:
    self._cutoff = float((value * self._cutoff_units).rescale(eV))
Clj.ewald_cutoff = property(_get_gcutoff, _set_gcutoff)
