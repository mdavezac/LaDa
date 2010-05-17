""" Classes to handle ESCAN wavefunctions. """


from ..opt.decorators import broadcast_result, make_cached
from ._extract import Extract

def rWavefunction(Object):
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
      result += dot( conjugate(self.down), multiply(operator, self.down) if operator != None else self.down)
    return result

  def braket(operator, ket):
    """ Returns <this wavefunction|operator|ket wavefuntion>.
    
        The operator is applied to the ket. If operator is None, then simply
        perform scalar product.
    """
    from numpy import conjugate, dot, multiply
    result = dot(conjugate(self.up), multiply(operator, ket.up) if operator != None else ket.up)
    if down != None: 
      result += dot(conjugate(self.down), multiply(operator, ket.down) if operator != None else ket.down)
    return result

def Wavefunction(rWavefunction):
  is_gspace = True
  """ True since this is a g-space wavefunction. """

  def __init__(self, index, eigenvalue, up, down = None, attenuation=None):
    """ Initializes a spinor. """
    super(Wavefunction, self).__init__(index, eigenvalue, up, down)
    self.attenuation = attenuation

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
      ket_ = multiply(ket.down, self.attenuation) i
      result += dot(conjugate(bra), multiply(operator, ket_) if operator != None else ket_)
    return result


class Wavefunctions(Extract):
  """ G-space wavefunctions. """
  is_gspace = True
  """ True if the wavefunctions are in gspace, False if  in realspace. """
  def __init__(self, *args, **kwargs)
    super(Wavefunctions, self).__init__(self, *args, **kwargs)


  def __getitem__(self, index): return self._wavefunction.__getitem(index)
  def __setitem__(self, index): raise RuntimeError("Wavefunctions cannot be modified.")
  def __len__(self): return self._wavefunctions.__len__(self)

  @property
  @make_cached
  def _wavefunctions(self):
    """ Creates list of Wavefuntion objects.

        This is makes instantiations "lazy", in that wavefunctions are
        resolved, eg read from disk, only when truly needed, not a
        instantiation per say.
    """
    result = []
    if self.is_spinor:
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(i, eig, self.raw_wfns[:,i,0],\
                                        self.raw_wfns[:,i,1], self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(i, eig, self.raw_wfns[self.inverse_indices,i,0],\
                                        self.raw_wfns[self.inverse_indices,i,1], self.attenuation) )
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(i, eig, self.raw_wfns[:,i,0],\
                                      self.raw_wfns[:,i,1], self.attenuation) )
    else: # no spin polarization.
      if self.is_krammer:
        for i, eig in enumerate(self.eigenvalues):
          if i % 2 == 0: # normal
            result.append( Wavefunction(i, eig, self.raw_wfns[:,i,0],\
                                        None, self.attenuation) )
          else:  # inverted
            result.append( Wavefunction(i, eig, self.raw_wfns[self.inverse_indices,i,0], \
                                        None, self.attenuation) )
          result.append(result[-1])
      else: # no krammer degeneracy
        for i, eig in enumerate(self.eigenvalues):
          result.append( Wavefunction(i, eig, self.raw_wfns[:,i,0],None, self.attenuation) )
          result.append(result[-1])

  @property
  @make_cached
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
  def raw_wfns(self):
    """ Raw wavefunction data. """
    return self._raw[0]

  @property
  def gvectors(self):
    """ G-vector values. """
    return self._raw[1]

  @property
  def attenuation(self):
    """ G-vector attenuation values. """
    return self._raw[2]

  @property
  def inverse_indices(self):
    """ Indices to -G vectors. """
    return self._raw[3]

  @property
  @make_cached
  @broadcast_result(attr=True, which=0)
  def _wavefunction_path(self): return self.solo().escan.WAVECAR

  @property
  def is_spinor(self):
    """ True if wavefunction is a spinor. """
    from . import soH
    return self.escan.potential == soH

  @property
  def is_krammer(self):
    """ True if wavefunction is a spinor. """
    from numpy.linalg import norm
    from . import soH
    return norm(self.escan.kpoint) < 1e-12

class rWavefunctions(Wavefunctions, rawdata):
  """ R-space wavefunctions. """
  is_gspace = False
  """ False since this is r-space wavefunctions. """
  def __init__(self, *args, **kwargs)
    super(rWavefunctions, self).__init__(self, *args, **kwargs)
    self._data = rawdata[0]
    """ Contains the raw data for r-space wavefunctions. """
    self.rvectors = rawdata[1]
    """ R-space vectors. """

