""" Defines Local Density of States.

    Both a functional (`Functional`) and a function (`ldos`) are defined which
    compute the Local Density of State.  The LDOS is obtaine as a Fourier
    transform of the wavefunctions to the input real-space positions. In this
    way, positions are not limited to the FFT mesh.

    The functional returns an extraction object which caches the LDOS to file.
"""
__docformat__ = "restructuredtext en"
__all__ = ['ldos', 'Extract', 'Functional']

from .kescan import KEscan, Extract as KExtract
from ..opt import make_cached, FileCache

class _ldosfunc(object):
  """ Local density of states for a given set of positions within a given structure. """
  def __init__(self, eigenvalues, rs):
    """ Initializes a local density of state functor. 
    
        :Parameters:
          eigenvalues : numpy array
            Vector of eigenvalues. Can be signed with a unit, or without, in
            which case eV is assumed.
          rs : numpy array
            Matrix of densities per real-space position(row) and per band(column). 
    """
    from numpy import sqrt, pi
    from quantities import eV
    self.eigenvalues = eigenvalues.copy()
    """ Vector of eigenvalues. """
    self.rs = rs.copy()
    """ Matrix of densities per real-space position(row) and per band(column). """
    if not hasattr(self.eigenvalues, 'rescale'): self.eigenvalues *= eV 
    else: self.eigenvalues = self.eigenvalues.rescale(eV)
    self._inv_sqrt_pi = 1e0/sqrt(pi)
    """ Normalization constant 1e0/sqrt(|pi|). 

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """

  def __call__(self, energy, sigma=0.1):
    """ Calls smearing function over densities.
    
        :Parameters: 
          energy : float or scalar array
            Energy at which to compute local density of states. This can be real
            number, in which case it should be in eV, or a numpy scalar with a
            unit (from quantity).
          sigma : float or scalar array
            Width at half maximum of the gaussian smearing. This can be real
            number, in which case it should be in eV, or a numpy scalar with a
            unit (from quantity).
    """
    from numpy import dot, exp, array
    from quantities import eV
    if not hasattr(sigma, 'rescale'): sigma *= eV
    else: sigma = sigma.rescale(eV)
    if not hasattr(energy, 'rescale'): energy = array(energy) * eV
    else: energy = energy.rescale(eV)
    y = array([ [(E - e)/sigma for e in energy] for E in self.eigenvalues])
    return dot(self.rs, self._inv_sqrt_pi/sigma * exp(-y*y))


def ldos(extractor, positions, raw=False):
  """ Local density of states at given positions.
  
      :Parameters:
        extractor 
          Output from a `KEscan` or a `Functional`  calculation.
        positions : nx3 array
          Array with positions in real space where to compute LDOS.
        raw : boolean
          Whether to return the raw data or the LDOS itself, i.e. a function of
          the energy.
  """
  from numpy import tensordot, multiply, conjugate, exp, concatenate, array, rollaxis

  assert isinstance(extractor, KExtract),\
         ValueError('extractor argument should be KExtract isntance.')

  extractor = extractor.copy(unreduce=False)
  istr, ostr = extractor.input_structure, extractor.structure
  normalization = 0e0
  perpoint = []
  for i, equivs in enumerate(extractor.functional.kpoints.iter_equivalents(istr, ostr)):
    # computes all positions including symmetry equivalents.
    # Since we expect fewer real space points than fourrier space points,
    # symmetric equivalents are added to real-space positions. 
    # See for instance "Electronic Structure", Richard M. Martin, first
    # edition, chapter 4 section 5.
    extract = extractor[i]
    equivs = [u for u in equivs]
    operators = [op.inverse for index, m, k, op in equivs]
    all_positions = array([op(u) for op in operators for u in positions])
    multiplicities = [m for index, m, k, op in equivs]
    normalization += sum(multiplicities)

    # creates array which may include krammer degenerate.
    if extract.is_krammer:
      inverse = conjugate(extract.raw_gwfns[extract.inverse_indices,:,:])
      gwfns = concatenate((extract.raw_gwfns, inverse), axis=1)
    else: gwfns = extract.raw_gwfns
    # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
    v = exp(-1j * tensordot(all_positions, extract.gvectors, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    rspace = tensordot(v, gwfns, ((1),(0)))
    rspace = multiply(rspace, conjugate(rspace)).real
    # Sum over spin channels if necessary.
    if not extract.is_spinor: rspace = rspace[:,:,0]
    else: rspace = rspace[:,:,0] + rspace[:,:,1]
    # Sum degenerate states if necessary.
    if extract.is_krammer:
      assert rspace.shape[1] % 2 == 0
      # sum krammer degenerate states together since same eigenvalue.
      rspace = rspace[:,:rspace.shape[1]//2] + rspace[:,rspace.shape[1]//2:]
    
    # sum over equivalent kpoints. 
    N = len(positions)
    if abs(multiplicities[0] - 1e0) > 1e-12: rspace[:N, :] *= multiplicities[0]
    for j, m in enumerate(multiplicities[1:]):
      if abs(m - 1e0) > 1e-12: rspace[:N, :] += m * rspace[(j+1)*N:(j+2)*N, :]
      else: rspace[:N, :] += rspace[(j+1)*N:(j+2)*N, :]

    # append to reduced kpoint ldos list.
    perpoint.append(rspace[:N,:].copy())

  # normalize results and concatenate.
  result = rollaxis(array(perpoint), 0,-1) / float(normalization)
  result = array([u.flatten() for u in result])
  
  return result if raw else _ldosfunc(extractor.eigenvalues, result)



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
    from ldos import ldos as outer_ldos
    return outer_ldos(self, self.positions, raw=True)

  @property
  def positions(self):
    """ Positions for which to compute LDOS. """
    from numpy import array
    if getattr(self.functional, 'positions', None) == None:
      return array([a.pos for a in self.structure.atoms])
    if not hasattr(self.functional.positions, '__call__'): return self.functional.positions
    return self.funtional.positions(self.structure)

  @property
  @make_cached
  def ldos(self):
    """ Local density of states for ``positions``. """
    return _ldosfunc(self.copy(unreduce=False).eigenvalues.flatten(), self.raw_ldos)
   
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
    out = super(Functional, self).__call__(*args, **kwargs)
    result = self.Extract(parent=out)
    result.ldos # computes and saves ldos.
    return result



      
