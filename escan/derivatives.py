""" Numerical energy derivatives. """

def reciprocal( escan, structure, direction, outdir = None, comm = None, order = 1, \
                nbpoints = None, stepsize = 1e-2, lstsq = None, **kwargs ):
  """ Computes effective mass for a given direction.

      @warning: escan.nbstates (or nbstates if passed as arg) must be the
        *spin-polarized* number of electrons, whatever escan thinks.
      @param structure: The structure for wich to compute effective masses.
      @type  structure: crystal.Structure
      @param escan: Emiprical pseudo-potential functional.
      @type  escan: Escan
      @param direction: direction for which to compute derivatives.
      @type  direction: 3-tuple 
      @param order: Highest order derivative to perform. Defaults to 1.
      @type  order: int
      @param nbpoints: Number of points with wich to compute derivatives.
        Should be at least order + 1. Default = order + 1. 
      @type  nbpoints: int
      @param stepsize: Distance between interpolation points. Default = 1e-2.
        Units of 2S{pi}/a, with C{a=structure.scale}.
      @type  stepsize: float
      @param kwargs: Extra parameters which are passed on first to escan (if escan
         object has an attribute with the same name), then to the linear least
         square fit method. Note that this will *not* change the external escan
         object.  This function is stateless. 
      @param lstsq: Linear least square method. The first two parameters should
         be same as numpy.linalg.lstsq. Other parameters can be passed as extra
         parameters. Defaults to numpy.linalg.lstsq.
      @return: Same as return from lstsq.
  """
  from os import getcwd
  from os.path import join
  from shutil import rmtree
  from copy import deepcopy
  from math import pow, pi, factorial
  from numpy import zeros, array, dot
  from numpy.linalg import norm, lstsq as np_lstsq
  from boost.mpi import world
  from ..physics import a0, Hartree

  # check input arguments
  order = int(order)
  if order <= 0: raise ValueError, "order should be at least one (%i)." % ( order )
  if nbpoints == None: nbpoints = order + 1
  else: nbpoints = int(nbpoints)
  if nbpoints < order + 1: raise ValueError, "The number of points shoule be at least order+1."
  direction = array(direction, dtype="float64")
  if norm(direction) < 1e-12: "Direction cannot be 0."
  direction /= norm(direction)
  if lstsq == None: lstsq = np_lstsq
  if comm == None: comm = world
  if outdir == None: outdir = getcwd()
  vffrun = kwargs.pop("vffrun", escan.vffrun)
  genpotrun = kwargs.pop("genpotrun", escan.genpotrun)

  # creates parameters matrix.
  parameters = zeros(shape=(nbpoints, order+1), dtype="float64") 
  parameters[:,0] = 1 # zero order terms are all 1.
  start = 0e0 
  if nbpoints % 2 == 0: start = 0.5 * stepsize
  units = 2e0 * pi * a0("A") / structure.scale 
  for i in range(0, (nbpoints-nbpoints%2) / 2):
    s = ( stepsize * float(i+1) + start ) * units
    parameters[2*i,   1:] = array([pow( s, n)/float(factorial(n)) for n in range(1, order+1)])
    parameters[2*i+1, 1:] = parameters[2*i,   1:]
    for j in range(1, order+1): 
      if j % 2 == 1: parameters[2*i+1, j] *= -1e0

  # first computes vff and genpot unless given.
  if genpotrun == None or vffrun == None: 
    vffout = escan( structure, outdir=outdir, do_escan=False, genpotrun=genpotrun,\
                    vffrun=vffrun, comm = comm, **kwargs )
    if genpotrun == None: genpotrun = vffout
    if vffrun == None: vffrun = vffout

  # now performs all calculations.
  measurements = None
  for i in range(0, nbpoints): 
    out = escan\
          (\
            structure, 
            outdir = join(outdir, "%i-s=%f" % (i, parameters[i,1] )),
            vffrun = vffrun, genpotrun = genpotrun, do_escan = True, 
            kpoint = escan.kpoint + direction * parameters[i, 1] / units,
            comm = comm,
            **kwargs
          )
    eigenvalues = out.eigenvalues.copy()
    eigenvalues.sort()
    if measurements == None: 
      measurements = zeros(shape=(nbpoints, eigenvalues.size), dtype="float64")
    measurements[i,:] = eigenvalues / Hartree("eV")
    
  # finally, performs least-square fit and returns evrything.
  result = lstsq( parameters, measurements )

  return result

