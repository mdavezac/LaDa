""" Classes to handle ESCAN wavefunctions. """


from ..opt.decorators import broadcast_result, make_cached
from ._extract import Extract

class Wavefunctions(Extract):
  def __init__(self, *args, **kwargs)
    super(Wavefunctions, self).__init__(self, *args, **kwargs)


  def __getitem__(self, index): pass
  def __setitem__(self, index): pass
  def __len__(self, index): pass

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _raw(self): pass

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _wavefunction_path(self): return self.solo().escan.WAVECAR
