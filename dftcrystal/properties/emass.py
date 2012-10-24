""" Computes effective masses for any structure. """
from ...tools.makeclass import makeclass, makefunc
from ...tools import make_cached
from .extract import Extract as ExtractBase
from . import Properties

class Extract(ExtractBase):
  """ Extraction for effective mass. """
  def __init__(self, *args, **kwargs):
     super(Extract, self).__init__(*args, **kwargs)

  @property 
  def success(self):
    """ True if successful run. 
    
        Checks this is an effective mass calculation.
    """
    try: self._details
    except: return False
    return ExtractBase.success.__get__(self)

  @property
  @make_cached
  def _details(self):
    """ Parameters when calling the effective mass routine. """
    from re import compile
    from ...misc import exec_input
    from ...error import GrepError
    start = compile(r'^#+ EMASS DETAILS #+$')
    end = compile(r'^#+ END EMASS DETAILS #')
    with self.__stdout__() as file:
      lines = None
      for line in file:
        if start.match(line): lines = ""; break
      if lines is None: raise GrepError('Could not find call parameters.')
      for line in file:
        if end.match(line):  break
        lines += line

    input = exec_input(lines)
    return { 'center': input.center,
             'polarpoints': input.polarpoints,
             'nbpoints': input.nbpoints,
             'range': input.range }
  @property
  def center(self): return self._details['center']
  @property
  def polarpoints(self): return self._details['polarpoints']
  @property
  def nbpoints(self): return self._details['nbpoints']
  @property
  def range(self): return self._details['range']
  @property
  def kpoints(self): return self.bandstructure.kpoints
  @property
  def eigenvalues(self): return self.bandstructure.eigenvalues

  @staticmethod
  def _orders(orders):
    """ Order up to which taylor coefficients should be computed. """
    result = orders 
    if result is None: result = [0, 2]
    if not hasattr(result, '__iter__'): result = [result]
    return result

  def parameters(self, orders):
    """ Parameters for the least-square-fit approach. """
    from numpy import zeros, pi, dot
    from numpy.linalg import inv
    from math import factorial

    # standardize order input. Now a list.
    order = self._orders(orders)
    # computes number of free parameters, accounting for analyticity.
    u = 0
    for order in orders:
      for p in xrange(order+1):
        for j in xrange(p+1): u += 1
    # creates parameter array
    result = zeros((len(self.kpoints), u), dtype="float64")
    # populates it.
    u = 0
    recipcell = inv(self.structure.cell).T
    kpoints = dot(recipcell, (self.kpoints - self.center).T).T * 2.0 * pi      \
              / self.structure.scale 
    stepsize = zeros(len(self.kpoints), dtype="float64")
    for order in orders:
      for p in xrange(order+1):
        for q in xrange(p+1):
          # constructs array of (x_j-x_0)^q_i, with multiplication over i
          # representing the cartesian coordinates and j the kpoints.
          stepsize[:]                                                          \
            = 1. / (factorial(order-p)*factorial(p-q)*factorial(q))
          for i, j in zip([order-p, p-q, q], xrange(3)): 
            if i != 0: stepsize *= kpoints[:,j]**i
          # fill in the parameters for this expansion.
          result[:, u] = stepsize 
          u += 1
    return result

  def fit(self, orders=None):
    """ Returns all fitted taylor coefficients. """
    from numpy import sum, argmin
    from numpy.linalg import lstsq
    if orders is None: orders = [0, 2]
    # computes measurements
    eigenvalues = self.eigenvalues
    if 0 not in orders: # recenter.
      index = argmin(sum((self.kpoints - self.center)**2, axis=1))
      if len(eigenvalues.shape) == 3: # spin  polarized
        eigenvalues[0] -= eigenvalues[0, index]
      else: # spin unpolarized
        eigenvalues -= eigenvalues[index]
    return lstsq(self.parameters(orders), eigenvalues.magnitude)

  def inverse_effective_mass(self, orders=None):
    """ Computes inverse effective mass tensors. """
    orders = self._orders(orders)
    startindex = 0
    if 0 in orders: startindex += 1
    if 1 in orders: startindex += 3
    if len(self.eigenvalues.shape) == 3

class _OnFinish(object):
  """ Called when effective mass calculation finishes. 
  
      Adds some data to the calculation so we can figure out what the arguments
      to the call.
  """
  def __init__(self, previous, outdir, center, polarpoints, nbpoints, range):
    self.previous = previous
    self.outdir = outdir
    self.center = [u for u in center]
    self.polarpoints = polarpoints
    self.nbpoints = nbpoints
    self.range = range
  def __call__(self, *args, **kwargs):
    from os.path import join
    # first calls previous onfinish.
    if self.previous is not None: self.previous(*args, **kwargs)
    # then adds data
    header = ''.join(['#']*20)
    with open(join(self.outdir, 'prop.out'), 'a') as file: 
      file.write('{0} {1} {0}\n'.format(header, 'EMASS DETAILS'))
      file.write('center={0!r}\n'.format(self.center))
      file.write('polarpoints={0!r}\n'.format(self.polarpoints))
      file.write('nbpoints={0!r}\n'.format(self.nbpoints))
      file.write('range={0!r}\n'.format(self.range))
      file.write('{0} END {1} {0}\n'.format(header, 'EMASS DETAILS'))

def iter_emass( functional, structure=None, outdir=None, workdir=None,
                center=None, polarpoints=5, nbpoints=5, range=0.1, **kwargs ):
  """ Computes effective mass for a given kpoint. 
  
      Computes Taylor expansion of a Hamiltonian around a given k-point.
      The expansion will be obtained by fitting an n-th order expansion to a
      set of eigenvalues in the neighberhood of the kpoint. In practice, a
      polar grid is created as a set of "band-structure" paths between
      diametrically opposed end-points. The eigenvalues are obtained for each
      path using CRYSTAL_'s BAND keyword (in properties program). 

      :param functional:
          The effective masses for the given k-point are computed for a given
          fixed Hamiltonian. This fixed hamiltonian must first be generated
          from a self-consistent electronic minization calculation. For
          convenience, the fixed Hamiltonian can be either be given  on input
          or  computed within this function. As such, ``functional`` can be one
          of two types.

          - :py:class:`~lada.dftcrystal.functional.Functional` object: a
             self-consistent calculation is first performed. The ``structure``
             argument cannot be None.
          - extraction object: e.g. as returned by a call to a
            :py:class:`~lada.dftcrystal.functional.Functional` object. It
            should point to the self-consistent calculation on which the
            effective mass calculations will be based. The ``structure``
            argument is ignored.

      :param structure:
        If the argument ``functional`` is indeed a CRYSTAL functional, then
        this routine will first perform as self-consistent electronic
        minization (or whatever the functional does) in order to obtain a
        Hamiltonian. This argument is the structure for which to compute the
        Hamiltonian. If the argument ``functional`` is an extraction object,
        than ``structure`` is ignored.

      :param outdir: Output directory.
      :param workdir: Directory where calculations are performed.
      :param center: k-point for which to perform taylor expansion.
      :type  center: 3d-vector in fractional coordinates
      :param int polarpoints:
        Number of points on a half-circle. These points and their diametrical
        oppposite define the "band-structure" lines for which eigenvalues are
        computed. 
      :param int nbpoints: 
        Number of points per "band-structure" lines.
      :param float range:
        
  """
  from numpy import mgrid, concatenate, pi, sin, cos, outer, ones, array, zeros
  from ...error import input as InputError
  from . import Properties

  # default kpoint is gamma.
  if center is None: center = zeros(3, dtype='float64')

  # If has an 'iter' function, then calls it. 
  if hasattr(functional, 'iter'): 
    if structure is None:
      raise InputError( 'If the first argument to iter_emass is a functional, '\
                        'then a structure must also be given on which to '     \
                        'apply the CRYSTAL functional.' )
    for input in functional.iter( structure, outdir=outdir, workdir=workdir,
                                  **kwargs ):
      if getattr(input, 'success', False): continue
      elif hasattr(input, 'success'):
        yield Extract(outdir)
        return
      yield input 
  # if is callable, then calls it.
  elif hasattr(functional, '__call__'):
    input = functional(structure, outdir=outdir, workdir=workdir, **kwargs)
  # otherwise, assume it is an extraction object.
  else: input = functional

  # check that self-consistent run was successful.
  if not input.success:
    yield Extract(outdir, input=input)
    return

  # creates Properties functional
  functional = Properties(input=input)
  # creates polar grid
  theta = mgrid[0:pi:pi/float(polarpoints)]
  cos_theta = cos(theta)
  sin_theta = sin(theta)
  x0 = range * outer(cos_theta, sin_theta).flatten()
  y0 = range * outer(sin_theta, sin_theta).flatten()
  z0 = range * outer(ones(len(cos_theta)), cos_theta).flatten()
  # The "grid" of lines is polar and centered at the kpoint.
  # In practice, it translates into a list of diametrically opposite points.
  # points: points on one side of ball centered at [0, 0, 0]
  points = concatenate((x0[:, None], y0[:, None], z0[:, None]), axis=1)
  # points contains replicas of the same point at the pole. we know make sure
  # each point is unique. we also make sure that diametrically opposed points
  # are avoided (there shouldn't be any, but whatever).
  unique = [points[0]]
  for point in points[1:]: 
    if any(sum((point-p)**2) < 1e-8 for p in unique): continue
    if any(sum((point+p)**2) < 1e-8 for p in unique): continue
    unique.append(point)
  unique = array(unique)
  # endpoints: pairs of diametically opposite points on a ball centered at kpoint.
  endpoints = concatenate(( (center+unique)[:,None, :],
                            (center-unique)[:, None, :]), axis=1)
  functional.band = endpoints
  # number of points per line should be nbpoints.
  functional.band.nsub = len(functional.band) * nbpoints

  # now call functional.
  for u in functional.iter(outdir=outdir, workdir=workdir):
    if getattr(u, 'success', False): continue
    if hasattr(u, 'success'): yield u; return
    # modify onfinish so that call arguments are added to the output file.
    onfinish = _OnFinish( u.onfinish, outdir, center, polarpoints, nbpoints,
                          range)
    u.onfinish = onfinish
    yield u
    
  yield Extract(outdir, input=input)

iter_emass.Extract = Extract
""" Extractor class for effective mass. """
EMass = makeclass( 'EMass', Properties, iter_emass, None,
                    module='lada.dftcrystal.properties.emass',
                    doc='Functional form of the '                              \
                        ':py:class:`lada.dftcrystal.properties.emass'          \
                        '.iter_emass` method.' )

# Function call to effective mass. No iterations. returns when calculations are
# done or fail.
effective_mass = makefunc( 'effective_mass', iter_emass,
                           'lada.dftcrystal.properties.emass' )
