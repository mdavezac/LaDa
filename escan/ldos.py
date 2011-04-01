""" Defines Local Density of States. """
__docformat__ = "restructuredtext en"
__all__ = ['ldos', 'Extract', 'Functional']

from . import KEscan, KExtract
from ..opt import make_cached, FileCache


class _ldosfunc(object):
  """ Local density of states for a given set of positions within a given structure. """
  def __init__(self, eigs, rs):
    from numpy import multiply
    self.eigs, self.rs = -multiply(eigs, eigs), rs.copy()
  def __call__(self, e, sigma=0.1):
    from numpy import dot, pi, exp
    return dot(self.rs, 1e0/sqrt(pi)/sigma * exp(-self.eigs/sigma/sigma))


def ldos(extract, positions, raw=False):
  """ Local density of states from previous calculation """
  from numpy import zeros, tensordot, multiply, conjugate, exp

  extractors = extract if hasattr(extract, "__iter__") else [extract]

  result = zeros(positions.shape[0], dtype="float64")

  for extract in  extractors:
    # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
    v = exp(-1j * tensordot(positions, extract.gvectors, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    rspace = tensordot(v, extract.raw_gwfns, ((1),(0)))
    # reduce across processes
    rspace = extract.comm.reduce(rspace, lambda x,y: x+y)
  
    if extract.is_krammer:
      rspace2 = rspace 
      rspace = zeros( (rspace.shape[0], rspace.shape[1]*2, rspace.shape[2]), dtype="complex64")
      rspace[:,::2,:] = rspace2
      cj = extract.raw_gwfns[extract.inverse_indices,:,:].conjugate()
      rspace2 = tensordot(v, cj, ((1),(0)))
      rspace2 = extract.comm.reduce(rspace, lambda x,y: x+y)
      rspace[:,1::2,:] = rspace2
    rspace = multiply(rspace, conjugate(rspace))
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
      assert len(kwargs) == 0 and len(args) == 0, \
             ValueError('Use of parent is exclusive')
    KExtract.__init__(self, *args, **kwargs)
    if parent != None: self.__dict__.update(parent.__dict__)
  
  @property
  @FileCache('LDOSCAR')
  def raw_ldos(self):
    """ Raw Local density of states for given sets of positions. """
    from .ldos import ldos as outer_ldos
    return outer_ldos(self, self.positions, raw=True)

  @property
  @make_cached
  def ldos(self):
    """ Local density of states for `positions`. """
    return _ldosfunc(self.eigenvalues, self.raw_ldos)
   
  @property
  def positions(self):
    """ Positions for which to compute LDOS. """
    from numpy import array
    if self.functional.positions == None: return array([a.pos for a in self.structure.atoms])
    if not hasattr(self.functional.positions, '__call__'): return self.functional.positions
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
    super(Functional, self).__init__(**kwargs)
    self.positions = kwargs.pop('positions', None)
    """ Callable returning positions for local density of states.

        Should be None (in which case atomic positions are used) or a
        pickleable callable which takes a structure and returns the positions
        for which to compute LDOS.
    """

  def __call__(self, *args, **kwargs):
    """ Calls KEscan, and then calls LDOS itself. 

        All parameters are passed on to escan.
    """ 
    print "BEFORE call"
    out = super(Functional, self).__call__(*args, **kwargs)
    print "BEFORE extract"
    result = self.Extract(parent=out)
    print "BEFORE LDOS"
    result.ldos # computes and saves ldos.
    print "AFTER LDOS"
    return result



      
