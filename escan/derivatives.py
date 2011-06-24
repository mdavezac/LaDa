""" Numerical energy derivatives. """
__docformat__ = "restructuredtext en"
from .kpoints import KPoints, _reduced_grids_factory

class DDPoints(KPoints):
  """ Points to compute Derivate of the Dispersion. """
  def __init__(self, direction, center = None, order = 2, nbpoints = 0, stepsize=1e-2, relax=True):
    """ Initializes the dispersion derivative k-point object. 


        :Parameters:
          direction : 3-tuple 
            Direction for which to compute derivatives.
          center : 3-tuple
            Point at which to take derivative.
          order : int
            Order of the derivative.
          nbpoints : int
            Number of points to use in computing derivative.
          stepsize : float
            Distance between interpolation points. Default = 1e-2.
            Units of ``2|pi|/a``, with ``a=structure.scale``.
          relax : bool
            Whether kpoint should be relaxed from input to output vff
            structure. True by default.

        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.direction = direction
    """ Direction for which to compute derivatives. """
    self.center = center
    """ Point at which to take derivative. """
    self.order = order
    """ Order of the derivative. """
    self.nbpoints = nbpoints
    """ Number of points to use in computing derivative. """
    self.stepsize = stepsize
    """ Distance between interpolation points. 

        Units of ``2|pi|/a``, with ``a=structure.scale``.
        
        .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI
    """
    self.relax = relax
    """ Whether to deform kpoints to the relaxed structure. """

  @property
  def parameters(self): 
    """ Parameters for derivation. """
    assert hasattr(self, "_parameters"),\
           RuntimeError("parameters attribute cannot be accessed befor calculation.")
    return self._parameters

  def _mnk(self, input, output):
    """ Yields lines of k-points to perform numerical derivation. """
    from math import factorial
    from numpy import array, dot, zeros, pi
    from numpy.linalg import norm, inv
    from quantities import angstrom
    from ..physics import a0

    assert norm(self.direction) > 1e-12, ValueError("Direction cannot be null.")
    assert self.order > 0, ValueError("Derivation order cannot be zero.")
    assert self.stepsize > 0, ValueError("Stepsize must be positive.")
    assert abs(output.scale) > 1e-12, ValueError("scale is zero in structure.")
    assert self.order > 0, ValueError("Order of derivative should be positive.")

    nbpoints  = max(self.order+1, self.nbpoints)
    direction = array(self.direction, dtype="float64") 
    center    = array(self.center, dtype="float64") if self.center != None\
                else zeros((3,), dtype="float64")
    if self.relax:
      deformation = dot(inv(output.cell.T), input.cell.T)
      direction = dot(deformation, direction)
      center = dot(deformation, center)
    direction /= norm(self.direction) # normalizes direction.

    # yields vector at point derivation for odd numbers of points.
    if nbpoints % 2 == 1: yield 1, center

    # yields all other points.
    start = 1 if nbpoints % 2 == 1 else 0.5
    parameters = zeros(shape=(nbpoints, self.order+1), dtype="float64") 
    parameters[:,0] = 1
    units = 2e0 * pi * a0.rescale(angstrom) / output.scale 
    for i in range(0, nbpoints - nbpoints%2): 
      # computes position on derivation line.
      s = (1 if i % 2 == 0 else -1) * self.stepsize * (i//2+start)
      # yields reciprocal space vector where to do calculation.
      yield 1, center + direction * s
      # sets up regression parameters.
      parameters[i + nbpoints%2, 1:] = [ pow(s*units, n)/float(factorial(n))\
                                         for n in range(1, self.order+1) ]

    # saves parameters for later use.
    self._parameters = parameters

  def __repr__(self):
    result = '{0.__class__.__name__}({1}'.format(self, repr(self.direction))
    do_key = self.center == None
    if self.center != None: result += ', {0}'.format(repr(self.center))
    if self.order != 2:
      result += ', {1}{0}'.format(repr(self.relax), 'relax=' if do_key else '') 
    else: do_key = True
    if self.nbpoints != 0:
      result += ', {1}{0.nbpoints}'.format(self, 'nbpoints=' if do_key else '') 
    else: do_key = True
    if self.stepsize != 1e-2:
      result += ', {1}{0.stepsize}'.format(self, 'stepsize=' if do_key else '') 
    else: do_key = True
    return result + ')'

ReducedDDPoints    = _reduced_grids_factory('ReducedDDPoints', DDPoints)


def reciprocal( escan, structure, outdir = None, comm = None, direction=(0,0,1), order = 1, \
                nbpoints = None, stepsize = 1e-2, center = None, lstsq = None, **kwargs ):
  """ Computes effective mass for a given direction.

      :Parameters:
        escan : `Escan` or `KEscan`
          Emiprical pseudo-potential functional.
        structure : `crystal.Structure`
          The structure for wich to compute effective masses.
        outdir : str
          Directory where to save results of calculation.
        comm : `lada.mpi.Communicator` or None
          MPI communicator containing processes with which to perform
          calculation.
        direction : 3-tuple 
          direction for which to compute derivatives.
        order : int
          Highest order derivative to perform. Defaults to 1.
        nbpoints : int
          Number of points with wich to compute derivatives.
          Should be at least order + 1. Default = order + 1. 
        stepsize : float
          Distance between interpolation points. Default = 1e-2.
          Units of ``2|pi|/a``, with ``a=structure.scale``.
        center : 3-tuple
          k-point where to take derivative. Units of ``2|pi|/a``, with
          ``a=structure.scale``.
        lstsq 
          Linear least square method. The first two parameters should
          be same as numpy.linalg.lstsq. Other parameters can be passed as extra
          parameters. Defaults to numpy.linalg.lstsq.
        kwargs 
          Extra parameters which are passed on first to escan (if escan
          object has an attribute with the same name), then to the linear least
          square fit method. Note that this will *not* change the external escan
          object.  This function is stateless. 

      :return: Same as return from lstsq.

      :warning: escan.nbstates (or nbstates if passed as arg) must be the
                *spin-polarized* number of electrons, whatever escan thinks.

      .. |pi|  unicode:: U+003C0 .. GREEK SMALL LETTER PI

  """
  from numpy import array, sort
  from numpy.linalg import lstsq as np_lstsq
  from quantities import hartree
  from .kescan import KEscan

  # takes care of default parameters.
  if not isinstance(escan, KEscan): escan = KEscan(escan=escan)
  if center == None: center = kwargs.pop("kpoint", escan.kpoint)
  center = array(center, dtype="float64")
  relax = kwargs.pop("do_relax_kpoint", escan.do_relax_kpoint)
  if outdir == None: outdir = "reciprocal"
  if lstsq == None: lstsq = np_lstsq

  # creates kpoints object.
  kpoints = ReducedDDPoints(direction, center, order, nbpoints, stepsize, relax)

  # performs calculations.
  out = escan(structure, outdir=outdir, comm=comm, kpoints=kpoints, **kwargs)
  # makes sure we have parameters. 
  for k in kpoints(out.input_structure, out.structure): continue

  # sorts eigenvalues at each kpoint and rescales to hartree.
  measurements = sort(out.eigenvalues.rescale(hartree), axis=1) 
  
  # finally, performs least-square fit and returns everything.
  result = lstsq( kpoints.parameters, measurements )

  return result
