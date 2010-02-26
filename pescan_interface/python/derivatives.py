""" Numerical energy derivatives. """

def reciprocal( structure, escan, vff, direction, order = 1, \
                nbpoints = None, stepsize = 1e-2, lstsq = None, **kwargs ):
  """ Computes effective mass for a given direction.

      @warning: escan.nbstates (or nbstates if passed as arg) must be the
        *spin-polarized* number of electrons, whatever escan thinks.
      @param structure: The structure for wich to compute effective masses.
      @type  structure: crystal.Structure
      @param escan: Emiprical pseudo-potential functional.
      @type  escan: Escan
      @param vff: Valence Force Field method.
      @type  vff: vff.Vff or vff.LayeredVff
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
  from os.path import join
  from shutil import rmtree
  from copy import deepcopy
  from math import pow, pi, factorial
  from numpy import zeros, array, dot
  from numpy.linalg import norm, lstsq as np_lstsq
  from ._escan import method, potential
  from ..physics import a0, Hartree

  # some sanity checks.
  order = int(order)
  if order <= 0: raise ValueError, "order should be at least one (%i)." % ( order )
  if nbpoints == None: nbpoints = order + 1
  else: nbpoints = int(nbpoints)
  if nbpoints < order + 1: raise ValueError, "The number of points shoule be at least order+1."
  if norm(direction) < 1e-12: "Direction cannot be 0."
  direction = array(direction, dtype="float64")
  direction /= norm(direction)
  if lstsq == None: lstsq = np_lstsq

  # saves mpicomm if necessary.
  mpicomm = escan.mpicomm
  if "mpicomm" in kwargs:
    mpicomm = kwargs["mpicomm"] 
    del kwargs["mpicomm"]
  # now copies escan
  escan = deepcopy(escan)
  # resets mpicomm 
  escan.mpicomm = mpicomm
  # sets to folded spectra
  escan.method = method.folded

  # sets other parameters.
  popthese = []
  for key in kwargs:
    if not hasattr(escan, key): continue
    setattr(escan, key, kwargs[key])
    popthese.append(key)
  for key in popthese: del kwargs[key]
  directory = escan.directory

  # number of states for which to perform calculations.
  nbstates = escan.nbstates

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


  # now performs all calculations.
  measurements = zeros(shape=(nbpoints, nbstates), dtype="float64")
  kpoint = escan.kpoint.copy()

  for i in range(0, nbpoints): 
    # computes kpoint
    escan.kpoint = kpoint + direction * parameters[i, 1] / units
    # checks for double/krammer mad degeneracy touble
    double_trouble = escan.potential != potential.spinorbit or norm(escan.kpoint) < 1e-12
    # in which case only half the eigenvalues are computed.
    if double_trouble: escan.nbstates = nbstates / 2
    # sets directory.
    escan.directory = join(join(directory, "recip_deriv"), "%i-s=%f" % (i, parameters[i,1] ))
    # actually computes stuff.
    result = escan(vff, structure)
    result.sort() # sorted eigs!
    if double_trouble: # in case escan tries to screw us up again.
      result = array([u for j in range(2) for u in result])
      escan.nbstates = nbstates
    # finally stores these results.
    measurements[i,:] = result / Hartree("eV")

  # finally, performs least-square fit and returns evrything.
  result = lstsq( parameters, measurements )

  if escan.destroy_directory == True: rmtree( join(directory, "recip_deriv") )

  return result

