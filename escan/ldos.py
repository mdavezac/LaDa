""" Defines Local Density of States. """
__docformat__ = "restructuredtext en"
__all__ = ['ldos', 'Extract', 'Functional']

from .kescan import KEscan, KExtract
from ..opt import make_cached, FileCached


class _ldosfunc(object):
  """ Local density of states for a given set of positions within a given structure. """
  def __init__(self, eigs, rs):
    self.eigs, self.rs = -np.multiply(eigs, eigs), rs.copy()
  def __call__(self, e, sigma=0.1):
    return np.dot(self.rs, 1e0/sqrt(np.pi)/sigma * np.exp(-self.eigs/sigma/sigma))


def ldos(extract, positions, raw=False):
  """ Local density of states from previous calculation """
  import numpy as np
  from .. import lada_with_mpi

  extractors = extract if not hasattr(extract, "__iter__") else [extract]

  result = np.zeros(positions.shape[0], dtype="float64")

  for extract in  extractors:
      # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
      v = np.exp(-1j * np.tensordot(positions, extract.gvectors, ((1),(1))))
      # computes fourrier transform for all wavefunctions simultaneously.
      rspace = np.tensordot(v, extract.raw_wfns, ((1),(0)))
      # reduce across processes
      if lada_with_mpi:
        from boost.mpi import reduce
        rspace = comm.reduce(rspace, lambda x,y: x+y)
  
      if extract.is_krammer:
        rspace2 = rspace 
        rspace = zeros( (rspace.shape[0], rspace.shape[1]*2, rspace.shape[2]), dtype="complex64")
        rspace[:,::2,:] = rspace2
        cj = extract.raw_wfns[extract.inverse_indices,:,:].conjugate()
        rspace2 = np.tensordot(v, cj, ((1),(0)))
        rspace2 = comm.reduce(rspace, lambda x,y: x+y)
        rspace[:,1::2,:] = rspace2
      rspace = np.multiply(rspace, np.conjugate(rspace))
      if not extract.is_spinor: rspace = rspace[:,:,0]
      else: rspace = rspace[:,:,0] + rspace[:,:,1]
      result += rspace
  result /= float(len(extractors))
  return result if raw else _ldosfunc(result)


class Extract(KExtract):
  """ Extraction routine for LDOS. """
  def __init__(self, *args, **kwargs): 
    """ Creates Extraction object. 


        All parameters are passed on to KExtract.__init__, unless
        ``parent`` is present. In that case, ``parent`` should be a KExtract
        object which will be copied. This way, we can add ldos specific
        properties to a KExtract object.
    """
    parent = kwargs.pop('parent', None)
    if parent != None:
      assert len(kwargs) == 0 and len(args) == 0, 
             ValueError('Use of parent is exclusive')
    KExtract.__init__(self, *args, *kwargs)
    if parent != None: self.__dict__.update(parent)
  
  @property
  @FileCached('LDOSCAR')
  def raw_ldos(self):
    """ Raw Local density of states for given sets of positions. """
    from . import ldos as outer_ldos
    return outer_ldos(self, self.positions, raw=True)

  @property
  @make_cached
  def ldos(self):
    """ Local density of states for `positions`. """
    return _ldosesfunc(self.eigenvalues, self.raw_ldos)
   
  @property
  def positions(self):
    """ Positions for which to compute LDOS. """
    return self.funtional.positions(self.structure)
      
  
  def iterfiles(self, **kwargs):
    """ Iterates through exportable files. 

        All parameters passed on to KExtract. 
        Adds LDOSCAR to export files.
    """ 
    from os.path import exists, join
    path = join(self.directory, 'LDOSCAR')
    if exists(path): yield path
    for file in KExtract.iterfiles(self, **kwargs): yield file

  @property 
  def success(self):
    """ True if successful run. """
    from os.path import join, exists
    if not exists(join(self.directory, 'LDOSCAR')): return False
    return KExtract.success.__get__(self)

class Functional(KEscan): 
  """ Functional to compute local density of states. """
  Extract = Extract
  """ Extraction object for LDOS. """
  def __init__(self, **kwargs):
    """ Initializes an LDOS functional. 
    
        :param kwargs: Any keyword argument that works for `KEscan`.
        :kwarg positions: callable which takes a structure and returns an
          array of positions where to perform ldos. Can be None, in which case,
          it defaults to the atomic positions. Must also be pickleable.
    """
    from quantities import angstrom
    self.positions = kwargs.pop(positions, None)
    """ Callable returning positions for local density of states.

        Should be None (in which case atomic positions are used) or a
        pickleable callable which takes a structure and returns the positions
        for which to compute LDOS.
    """

  def __call__(self, *args, **kwargs):
    """ Calls KEscan, and then calls LDOS itself. 

        All parameters are passed on to escan.
    """ 
    out = KEscan.__call__(self, *args, **kwargs)
    result = self.Extract(parent=out)
    result.ldos # computes and saves ldos.
    return result



      
