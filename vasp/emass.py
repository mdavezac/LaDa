""" Methods to compute effective masses and other derivatives. """
__docformat__ = "restructuredtext en"
__all__ = ["Extract", "iter_emass", 'EMass', 'effective_mass']
from .extract import Extract as ExtractDFT
from ..functools.makeclass import makeclass, makefunc
from ..functools import make_cached
from .functional import Vasp

class Extract(ExtractDFT):
  """ Extractor for reciprocal taylor expansion. """
  def __init__( self, outcar=None, input=None, lstsq=None,                     \
                order=None, rcond=-1, **kwargs ):
    """ Initializes the extraction object. 
    
        :Parameters:
          outcar : str or None
	    Path to OUTCAR file. Can also be the directory if the OUTCAR is
	    named "OUTCAR". Defaults to None, in which case it uses the
            current working directory.
          input : Extract-type
	    Extractor object for the calculation from which the charge density
            was obtained. Defaults to None, meaning whatever calculation is
            contained in the parent directory.
          lstsq : method
            Least-square-fit method from which to get taylor coefficients.
            Defaults to None, i.e. numpy's lstsq method.
          order : unsigned int or list 
            Order up to which to compute taylor coefficients.
    """
    from os.path import dirname, exists, join
    if exists(join(outcar, 'reciprocal')):
      outcar = join(outcar, 'reciprocal')
    super(Extract, self).__init__(outcar, **kwargs)
    self.input = input if input != None                                        \
                 else ExtractDFT(dirname(self.directory))
    """ Extraction object of charge density calculation. """
    self.lstsq = lstsq 
    """ Least-square-fit method from which to get Taylor coefficients. """
    self.order = order
    self.rcond = rcond

  @property
  def rcond(self): return self._rcond
  @rcond.setter
  def rcond(self, value):
    self._rcond = value
    if hasattr(self, "_parameters"): delattr(self, "_parameters")
    if hasattr(self, "_measurements"): delattr(self, "_measurements")
    if hasattr(self, "_fitresults"): delattr(self, "_fitresults")

  @property
  @make_cached
  def _details(self):
    """ Details of emass calculation. """
    extracted = ''
    with self.__outcar__() as file:
      for line in file: 
        if line == '################ EMASS DETAILS ################\n':
          break
      for line in file:
        if line == '################ END EMASS DETAILS ################':
          break
        extracted += '\n' + line 
    results = {}
    exec extracted in results
    return results

  @property
  def order(self):
    """ Order up to which taylor coefficients should be computed. """
    result = self._order 
    if result is None: result = self._details['order']
    if not hasattr(result, '__iter__'): result = [result]
    return result

  @order.setter
  def order(self, value):
    """ Order up to which taylor coefficients should be computed. """
    # case where we want to fit to specific orders only.
    if value is None: self._order = None
    elif hasattr(value, "__iter__"):
      for i in sorted(value):
        if i < 0: raise ValueError("Order cannot be less than 0.")
        try: nbpoints = self.nbpoints
        except: pass
        else:
          if i >= nbpoints:
            raise ValueError( "Size of kpoint mesh implies order cannot "      \
                              "be larger than {0}.".format(nbpoints-1) )
      self._order = list(sorted(value))
    # case where we want all orders up to given input.
    else: 
      if value < 0:
        raise ValueError("Order cannot be less than 0.")
      try: nbpoints = self.nbpoints
      except: pass
      else:
        if value >= nbpoints:
          raise ValueError( "Size of kpoint mesh implies order "               \
                            "cannot be larger than {0}.".format(nbpoints-1) )
      self._order = list(xrange(value+1))
    if hasattr(self, "_parameters"): delattr(self, "_parameters")
    if hasattr(self, "_measurements"): delattr(self, "_measurements")
    if hasattr(self, "_fitresults"): delattr(self, "_fitresults")

  @property 
  def lstsq(self): 
    """ Least square fitting method. """
    from numpy.linalg import lstsq as nplstsq
    if self._lstsq is None:  
      result = self._details['lstsq']
      return nplstsq if result is None else result
    return nplstsq if self._lstsq is None else self._lstsq
  @lstsq.setter
  def lstsq(self, value):
    self._lstsq = value
    if hasattr(self, "_fitresults"): delattr(self, "_fitresults")

  @property
  def nbpoints(self):
    """ Number of points per axis involved in calculation. """
    from math import pow
    return int(pow(len(self.kpoints), 1./3.))

  @property
  def center(self):
    """ Returns central kpoint. """
    return self.kpoints[0] if self.nbpoints % 2 == 0                           \
           else self.kpoints[len(self.kpoints)//2]

  @property
  def stepsize(self):
    """ Returns central kpoint. """
    N = len(self.kpoints)//2
    return (self.kpoints[N] - self.kpoints[N-1])[2]
    
  @property
  def parameters(self):
    """ Parameters for the least-square-fit approach. """
    from numpy import zeros, pi
    from math import factorial

    if not hasattr(self, "_parameters"): 
      # computes number of free parameters, accounting for analyticity.
      u = 0
      for order in self.order:
        for p in xrange(order+1):
          for j in xrange(p+1): u += 1
      # creates parameter array
      self._parameters = zeros((len(self.kpoints), u), dtype="float64")
      # populates it.
      u = 0
      kpoints = (self.kpoints - self.center) * 2.0 * pi                        \
                / self.structure.scale 
      stepsize = zeros(len(self.kpoints), dtype="float64")
      for order in self.order:
        for p in xrange(order+1):
          for q in xrange(p+1):
            # constructs array of (x_j-x_0)^q_i, with multiplication over i
            # representing the cartesian coordinates and j the kpoints.
            stepsize[:]                                                        \
              = 1. / (factorial(order-p)*factorial(p-q)*factorial(q))
            for i, j in zip([order-p, p-q, q], xrange(3)): 
              if i != 0: stepsize *= kpoints[:,j]**i
            # fill in the parameters for this expansion.
            self._parameters[:, u] = stepsize 
            u += 1
    return self._parameters

  @property
  def measurements(self):
    """ Measurements for the least-square fit approach. 
    
        There is currently no attempt at accounting for band-crossing.
        This is simply the eigenvalues reshaped so spin is unrolled.
    """
    if self.ispin == 1: return self.eigenvalues.magnitude
    shape = self.eigenvalues.shape[0] * self.eigenvalues.shape[1],             \
            self.eigenvalues.shape[2]
    return self.eigenvalues.reshape(*shape).magnitude

  @property 
  def fitresults(self):
    """ Returns all fitted taylor coefficients. """
    from numpy.linalg import lstsq
    if not hasattr(self, "_fitresults"): 
      lstsq = getattr(self, 'lstsq', lstsq)
      self._fitresults = lstsq(self.parameters, self.measurements, self.rcond)
    return self._fitresults
    
  @property
  def residues(self):
    """ Return residues for each band. """
    return self.fitresults[1]
  @property
  def tensors(self):
    """ Return tensors for each band. """
    from numpy import zeros
    from itertools import chain, permutations
    from quantities import angstrom
    all_xs = self.fitresults[0].T
    result = [[]]*(len(self.order))
    
    current_index, current_order = 0, 0
    # add results for order 0.
    if self.order[0] == 0:
      result[current_order] = all_xs[current_index].copy()                     \
                              * self.eigenvalues.units                         
      current_index += 1                                                       
      current_order += 1                                                       
    # add results for order 1.                                                 
    if 1 in self.order:                                                        
      result[current_order]                                                    \
           = [ i.copy() * self.eigenvalues.units * angstrom                    \
               for i in all_xs[current_index:current_index+3].T ]
                              
      current_index += 3
      current_order += 1

    # compute index ranges for subsequent orders.
    u, indices_range = current_index, [current_index]
    for order in self.order[current_order:]:
      for p in xrange(order+1):
        for q in xrange(p+1): u += 1
      indices_range.append(u)

    # loop over remaining orders.
    for order, startu in zip(self.order[current_order:], indices_range):
      # loop over results for each band.
      for band_x in all_xs:
        u = startu
        # create tensor from n vector.
        dummy = zeros([3]*order, dtype="float64")
        for p in xrange(order+1):
          for q in xrange(p+1):
            indices = [ [i]*j for i, j in
                         zip([0, 1, 2], [order-p, p-q, q]) if j != 0 ]
            indices = set(permutations(chain(*indices)))
            for index in indices:
              dummy[index] = band_x[u] if abs(band_x[u]) > 1e-12 else 0
            u += 1
        # add tensor to results for that order.
        dummy *= self.eigenvalues.units * angstrom**order
        result[current_order].append(dummy)
      current_order += 1

    # got it all.
    return result

  @property
  def inverse_mass_tensor(self):
    """ Returns inverse mass tensor, hopefully with right units. """
    from ..physics import electronic_mass, h_bar
    if 2 not in self.order:
      raise AttributeError( "Effective mass are not part "                     \
                            "of the current expansion." )
    result = self.tensors[self.order.index(2)]
    for i in xrange(len(result)): 
      result[i] = (result[i] / h_bar**2).rescale(1./electronic_mass)
    return result

  def iterfiles(self, **kwargs):
    """ Exports files from both calculations. """
    from itertools import chain
    for file in chain( super(Extract, self).iterfiles(**kwargs),
                       self.input.iterfiles(**kwargs) ):
      yield file
    
def iter_emass( vasp, structure, outdir=None, order=2, nbpoints=None,
               stepsize=1e-2, center=None, lstsq=None,
               emassparams=None, **kwargs ):
  """ Computes k-space taylor expansion of the eigenvalues up to given order.

      First runs a vasp calculation using the first input argument, regardless
      of whether a restart keyword argument is also passed. In practice,
      following the general LaDa philosophy of never overwritting previous
      calculations, this will not rerun a calculation if one exists in
      ``outdir``. 
      Second, a static non-self-consistent calculation is performed to compute
      the eigenvalues for all relevant kpoints.

      :param vasp: :py:class:`~vasp.Vasp` or derived functional.
      :type vasp: `vasp.Vasp`
      
      :param structure: The structure for wich to compute effective masses.
      :type structure: `~lada.crystal._cppwrapper.Structure`

      :param str outdir:
          Root directory where to save results of calculation. Calculations
	  will be stored in  "reciprocal" subdirectory of this input parameter.

      :param int order: Highest order derivative to perform. Defaults to 1.

      :param int nbpoints:
          Number of points (in a single direction) with wich to compute
          taylor expansion.  Should be at least order + 1. Default to order +
          1. Note that in the case of even nbpoints, an extra kpoint is added
          so that the center does get computed.

      :param float stepsize:
          Distance between interpolation points. Default = 1e-2.
          Units of ``2|pi|/a``, with ``a=structure.scale``.
      :param center:
          Central k-point of the taylor expansion. This should be given in
          **reciprocal** units (eg coefficients to the reciprocal lattice
          vectors). Default is None and means |Gamma|.
      :type center: 3 floats

      :param lstsq:
          Linear least square method. If None, defaults to numpy.linalg.lstsq_
      :type lstsq: callable

      :param dict emassparams: 
         Parameters for the (non-self-consistent) effective mass caclulation
         proper. For instance, could include pertubative spin-orbit
         (:py:attr:`~vasp.functional.Vasp.lsorbit`).

      :param kwargs:
          Extra parameters which are passed on to vasp, both for the initial
          calculation and the effective mass calculation proper.

      :return: Extraction object from which masses can be obtained.

      .. _numpy.linalg.lstsq::
        http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
      .. |pi|     unicode:: U+003C0 .. GREEK SMALL LETTER PI
      .. |Gamma|  unicode:: U+00393 .. GREEK CAPITAL LETTER GAMMA
  """
  from inspect import isfunction
  from copy import deepcopy
  from os import getcwd
  from os.path import join
  from numpy import array, dot, append
  from numpy.linalg import inv

  # save input details for printing later on.
  details = "order    = {0!r}\n"                                               \
            "nbpoints = {1!r}\n"                                               \
            "stepsize = {2!r}\n"                                               \
            "center   = {3!r}\n"                                               \
            .format( order, nbpoints, stepsize,                                \
                     None if center is None else list(center) )                
  if isfunction(lstsq):                                                        
    details =  "from {0.__module__} import {0.__name__} as lstsq\n"            \
               .format(lstsq) + details
  elif lstsq is None: 
    details += "lstsq    = None\n"
  else:
    details =  "from {0.__class__.__module__} import "                         \
                    "{0.__class__.__name__}\n".format(lstsq)                   \
               + details
    details += "lstsq    = {0!r}\n".format(lstsq)
 
  # takes care of default parameters.
  if center is None: center = kwargs.pop("kpoint", [0,0,0])
  center = array(center, dtype="float64")
  if outdir is None: outdir = getcwd()
  if nbpoints == None: nbpoints = order + 1
  if nbpoints < order + 1:
    raise ValueError( "Cannot compute taylor expansion of order {0} "          \
                      "with only {1} points per direction."                    \
                      .format(order, nbpoints) )

  # first runs vasp.
  for first in vasp.iter(structure, outdir=outdir, **kwargs): yield first
  if not first.success:
    raise RuntimeError('Input run, or relaxation, was not successfull.')

  # prepare second run.
  vasp = deepcopy(vasp)
  center = dot(inv(first.structure.cell).T, center)
  kpoints = [ (x, y, z) for x in xrange(nbpoints)                              \
                        for y in xrange(nbpoints)                              \
                        for z in xrange(nbpoints) ]                            
  vasp.kpoints = (array(kpoints, dtype="float64")                              \
                 - (nbpoints-1)/2.) * stepsize + center
  if nbpoints % 2 == 0: # adds central point to grid.
    vasp.kpoints = append(center[None,:], vasp.kpoints, axis=0)

  # and exectute it.
  kwargs = deepcopy(kwargs)
  kwargs['restart'] = first
  kwargs['nonscf']  = True
  kwargs['relaxation']  = None
  if emassparams is not None: kwargs.update(emassparams)
  directory = join(first.directory, "reciprocal")
  for u in vasp.iter(first.structure, outdir=directory, **kwargs):
    yield u

  # Add computational details here.
  with open(join(directory, 'OUTCAR'), 'a') as outcar:
    outcar.write('\n################ EMASS DETAILS ################\n')
    outcar.write(details)
    outcar.write('\n################ END EMASS DETAILS ################\n')

  yield iter_emass.Extract(first.directory)

iter_emass.Extract = Extract
""" Extractor class for the reciprocal method. """
EMass = makeclass( 'EMass', Vasp, iter_emass, None, module='lada.vasp.emass',
                    doc='Functional form of the '                              \
                        ':py:class:`lada.vasp.relax.iter_emass` method.' )

# Function call to effective mass. No iterations. returns when calculations are
# done or fail.
effective_mass = makefunc('effective_mass', iter_emass, 'lada.vasp.emass')
