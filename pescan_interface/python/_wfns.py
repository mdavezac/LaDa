""" Classes to handle ESCAN wavefunctions. """


from ..opt.decorators import broadcast_result, make_cached
from ._extract import Extract

class rWavefunction(object):
  is_gspace = False
  """ False since this is a r-space wavefunction. """

  def __init__(self, index, eigenvalue, up, down = None):
    """ Initializes a spinor. """
    self.index = index
    self.eigenvalue = eigenvalue
    self.up = up
    self.down = down

  def expectation_value(operator):
    """ Returns expectation value of this operator. 
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply
    result = dot( conjugate(self.up), multiply(operator, self.up) if operator != None else self.up)
    if down != None: 
      result += dot( conjugate(self.down), multiply(operator, self.down)\
                                           if operator != None else self.down)
    return result

  def braket(operator, ket):
    """ Returns <this wavefunction|operator|ket wavefuntion>.
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply
    result = dot(conjugate(self.up), multiply(operator, ket.up) if operator != None else ket.up)
    if down != None: 
      result += dot(conjugate(self.down), multiply(operator, ket.down)\
                                          if operator != None else ket.down)
    return result

class Wavefunction(rWavefunction):
  is_gspace = True
  """ True since this is a g-space wavefunction. """

  def __init__(self, *args, **kwargs):
    """ Initializes a spinor. """
    self.attenuation = kwargs["attenuation"] if "attenuation" in kwargs else None
    """ Attenuation coefficients of high-energy G-vectors. """
    if "attenuation" in kwargs: del kwargs["attenuation"]
    super(Wavefunction, self).__init__(*args, **kwargs)

  def expectation_value(operator, attenuate = True):
    """ Returns expectation value of this operator. 
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply
    if not attenuate: return super(Wavefunction, self).expectation_value(operator)
    u = multiply(self.up, self.attenuation) 
    result = dot( conjugate(u), multiply(operator, u) if operator != None else u)
    if down != None: 
      u = multiply(down, self.attenuation) 
      result += dot( conjugate(u), multiply(operator, u) if operator != None else u)
    return result

  def braket(operator, ket, attenuate = True):
    """ Returns <this wavefunction|operator|ket wavefuntion>.
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply
    if not attenuate: return super(Wavefunction, self).braket(operator, ket)
    bra = multiply(self.up, self.attenuation)
    ket_ = multiply(ket.up, self.attenuation)
    result = dot(conjugate(bra), multiply(operator, ket_) if operator != None else ket_)
    if down != None: 
      u = multiply(down, self.attenuation) 
      ket_ = multiply(ket.down, self.attenuation)  
      result += dot(conjugate(bra), multiply(operator, ket_) if operator != None else ket_)
    return result


def gtor_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  """ g-space to r-space fourrier transform of wavefunctions.
  
      @param wavefunctions: an numpy array of wavefunctions.
      @param rvectors: a two-dimensional array of r-space vectors, with each
        row a position. The return r-space wavefunctions will be given with
        respect to these points. Each process should have different r-space
        vectors. Otherwise use another implementation.
      param gvectors: a two-dimensional array of g-space vectors, with each
        row a (g-)position. The input wavefunctions should be given with
        respect to these points, in the same order, etc.
      @param comm: communicator over which the wavefunctions are distributed.
        The return wavefunctions will also be dirstributed over these
        processes.
      @params axis: axis over which the wavefunctions are deployed, eg axis
        of L{wavefunctions} which corresponds to L{gvectors}. Independent
        (though simultaneous, implementation wise) fourrier transform will be
        performed over this axis for all other axis. 0 by default.
      @type axis: integer
  """
  from sys import getrefcount
  import numpy as np
  from boost.mpi import broadcast, reduce

  result = None
  for node in range(comm.size):
    # sends rvectors from node to all
    r = broadcast(comm, rvectors, node)
    # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
    v = np.exp(-1j * np.tensordot(r, gvectors, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    dummy = np.tensordot(v, wavefunctions, ((1),(axis)))
    # reduce across processes
    if node == comm.rank: result = reduce(comm, dummy, lambda x,y: x+y, node)
    else: reduce(comm, dummy, lambda x,y: x+y, node)

  assert not np.any(np.isnan(result))
  return result


def rtog_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  """ r-space to g-space fourrier transform of wavefunctions.
  
      @param wavefunctions: an numpy array of wavefunctions.
      @param rvectors: a two-dimensional array of r-space vectors, with each
        row a position. The return r-space wavefunctions will be given with
        respect to these points. Each process should have different r-space
        vectors. Otherwise use another implementation.
      @param gvectors: a two-dimensional array of g-space vectors, with each
        row a (g-)position. The input wavefunctions should be given with
        respect to these points, in the same order, etc.
      @param comm: communicator over which the wavefunctions are distributed.
        The return wavefunctions will also be dirstributed over these
        processes.
      @params axis: axis over which the wavefunctions are deployed, eg axis
        of L{wavefunctions} which corresponds to L{gvectors}. Independent
        (though simultaneous, implementation wise) fourrier transform will be
        performed over this axis for all other axis. 0 by default.
      @type axis: integer
  """
  from sys import getrefcount
  import numpy as np
  from boost.mpi import broadcast, reduce, all_reduce

  assert not np.any(np.isnan(wavefunctions))
  assert not np.any(np.isnan(gvectors))
  assert not np.any(np.isnan(rvectors))
  result = None
  for node in range(comm.size):
    # sends rvectors from node to all
    g = broadcast(comm, gvectors, node)
    # computes all exponentials exp(-i r.g), with g in first dim, and r in second.
    v = np.exp(1j * np.tensordot(rvectors, g, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    # somehow, there is a problem with tensordot leading to nan numbers...
    assert not np.any(np.isnan(v))
    last = wavefunctions.ndim-1
    dummy = np.dot(wavefunctions.swapaxes(axis, last), v).swapaxes(last, axis)
    # reduce across processes
    assert not np.any(np.isnan(dummy))
    if node == comm.rank: result = reduce(comm, dummy, lambda x,y: x+y, node)
    else: reduce(comm, dummy, lambda x,y: x+y, node)


  assert not np.any(np.isnan(result))
  # gets normalization factor.
  norm = all_reduce(comm, rvectors.shape[0], lambda x,y:x+y)
  assert norm != 0
  return result/ float(norm)
