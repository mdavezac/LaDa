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
  def _raw(self):
    """ Reads and caches g-space wavefunction data. 
    
        This is a tuple described making up the return of
        L{read_wavefunctions<lada.escan._escan.read_wavefunctions>}
    """
    from ._escan import read_wavefunctions

    assert self.comm.size >= self.nnodes,\
           RuntimeError("Must read wavefunctions with at least "\
                        "as many nodes as they were written to disk.")
    if self.comm.size > self.nnodes:
      color = 0 if self.comm.rank < self.nnodes else 1
      local_comm = self.comm.split(color)
    else: color, local_comm = 0, self.comm
    if color == 1: return (,,,)
    return read_wavefunctions(self.escan, [i for i in range(self.nbstates)], local_comm)



  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _wavefunction_path(self): return self.solo().escan.WAVECAR
