""" Methods to compute effective masses and other derivatives. """
__docformat__ = "restructuredtext en"
__all__ = ["Extract", "reciprocal"]
from .extract import ExtractDFT

class Extract(ExtractDFT):
  """ Extractor for reciprocal taylor expansion. """
  def __init__(self, outcar=None, comm=None, input=None, **kwargs):
    """ Initializes the extraction object. 
    
        :Parameters:
          outcar : str or None
            Path to OUTCAR file. Can also be the directory if the OUTCAR is
            named "OUTCAR". Defaults to None, in which case it uses the current
            working directory.
          comm : None or `lada.mpi.Communicator`
            MPI Communicator grouping processes participating in the calculation.
          input : Extract-type
            Extractor object for the calculation from which the charge density
            was obtained. Defaults to None, meaning whatever calculation is
            contained in the parent directory.
          order : unsigned int or list 
            Order up to which to compute taylor coefficients.
    """
    from os.path import dirname
    from lada.vasp import Extract as VaspExtract
    super(Extract, self).__init__(outcar, comm=comm, **kwargs)
    self.input = input if input != None else VaspExtract(dirname(self.directory), comm=comm)
    """ Extractor object for the calculation from which the charge density was obtained. """

  def _order(self, value):
    """ Order up to which taylor coefficients should be computed. """
    # case where we want to fit to specific orders only.
    if hasattr(value, "__iter__"):
      for i in sorted(value):
        if i < 0: raise ValueError("Order cannot be less than 0.")
        try: nbpoints = self.nbpoints
        except: pass
        else:
          if i >= nbpoints:
            raise ValueError("Size of kpoint mesh implies order cannot be larger than {0}.".format(nbpoints-1))
      return list(sorted(value))
    # case where we want all orders up to given input.
    else: 
      if value < 0:
        raise ValueError("Order cannot be less than 0.")
      try: nbpoints = self.nbpoints
      except: pass
      else:
        if value >= nbpoints:
          raise ValueError("Size of kpoint mesh implies order cannot be larger than {0}.".format(nbpoints-1))
      return list(xrange(value+1))

  @property
  def nbpoints(self):
    """ Number of points per axis involved in calculation. """
    from math import pow
    return int(pow(len(self.kpoints), 1./3.)+1e-5)

  @property
  def center(self):
    """ Returns central kpoint. """
    return self.kpoints[0] if self.nbpoints % 2 == 0 else self.kpoints[len(self.kpoints)//2]
  @property
  def central_eigenvalues(self):
    """ Returns eigenvalues at center. """
    if self.ispin == 1:
      return self.eigenvalues[0] if self.nbpoints % 2 == 0 else self.eigenvalues[len(self.eigenvalues)//2]
    return self.eigenvalues[:,0,:] if self.nbpoints % 2 == 0 \
          else self.eigenvalues[:,len(self.eigenvalues)//2,:]

  def _parameters(self, order=2, kmag=True):
    """ Parameters for the least-square-fit approach. """
    from numpy import zeros, pi, abs
    from math import factorial

    # computes number of free parameters, accounting for analyticity.
    u = 0
    for o in self._order(order):
      for p in xrange(o+1):
        for j in xrange(p+1): u += 1
    # creates parameter array
    parameters = zeros((len(self.kpoints), u), dtype="float64")
    # populates it.
    u = 0
    kpoints = (self.kpoints - self.center) * 2.0 * pi / self.structure.scale 
    stepsize = zeros(len(self.kpoints), dtype="float64")
    for o in self._order(order):
      for p in xrange(o+1):
        for q in xrange(p+1):
          # constructs array of (x_j-x_0)^q_i, with multiplication over i
          # representing the cartesian coordinates and j the kpoints.
          stepsize[:] = 1. / (factorial(o-p)*factorial(p-q)*factorial(q))
          for i, j in zip([o-p, p-q, q], xrange(3)): 
            if i != 0: stepsize *= kpoints[:,j]**i
          # fill in the parameters for this expansion.
          parameters[:, u] = stepsize 
          u += 1
    return abs(parameters) if kmag else parameters

  def _measurements(self, order=2):
    """ Measurements for the least-square fit approach. 
    
        There is currently no attempt at accounting for band-crossing.
        This is simply the eigenvalues reshaped so spin is unrolled.
    """
    from numpy import rollaxis
    if 0 in self._order(order):
      if self.ispin == 1: return self.eigenvalues.magnitude
      shape = self.eigenvalues.shape[1], self.eigenvalues.shape[0]*self.eigenvalues.shape[2]
      return rollaxis(self.eigenvalues, 0, -1).reshape(*shape).magnitude
    if self.ispin == 1: return self.eigenvalues.magnitude - self.central_eigenvalues.magnitude
    shape = self.eigenvalues.shape[1], self.eigenvalues.shape[0]*self.eigenvalues.shape[2]
    return rollaxis(self.eigenvalues, 0, -1).reshape(*shape).magnitude\
           - self.central_eigenvalues.reshape(1,-1).magnitude

  def _fitresults(self, order=2, kmag=True, lstsq=None, **kwargs):
    """ Returns all fitted taylor coefficients. """
    from numpy.linalg import lstsq as np_lstsq
    if lstsq == None: lstsq = np_lstsq
    return lstsq(self._parameters(order, kmag), self._measurements(order), **kwargs)
    
  def residues(self, order=2, kmag=True, lstsq=None, tolerance=None, **kwargs):
    """ Return residues for each band. """
    return self._fitresults(order, kmag, lstsq, **kwargs)[1]

  def tensors(self, order=2, kmag=True, lstsq=None, tolerance=None, **kwargs):
    """ Return tensors listed by order. 
    
        Within each order, the tensors are listed as order then spin then
        bands if ispin == 2. In spin-degenerate cases, they are listed per band
        only.
    """
    from numpy import rollaxis
    fits = self._fitresults(order, kmag, lstsq, **kwargs)[0]
    if self.ispin == 1: return self._tensors_nospin(fits, order, tolerance)
    up, down = rollaxis(fits.reshape(self._fitresults.shape[0], 2, -1), 1)
    return list(zip(self._tensors_nospin(up, order, tolerance),
                    self._tensors_nospin(down, order, tolerance)))


  def _tensors_nospin(self, all_xs, order, tolerance):
    """ Return tensors for each band. """
    from numpy import zeros
    from itertools import chain, permutations
    from quantities import angstrom
    result = [[]]*(len(self._order(order)))
    
    current_index, current_order = 0, 0
    # add results for order 0.
    if self._order(order)[0] == 0:
      result[current_order] = all_xs[current_index].copy() * self.eigenvalues.units
      current_index += 1
      current_order += 1
    # add results for order 1.
    if 1 in self._order(order):
      result[current_order] = [ i.copy() * self.eigenvalues.units * angstrom\
                                for i in all_xs[current_index:current_index+3].T ]
                              
      current_index += 3
      current_order += 1

    # compute index ranges for subsequent orders.
    u, indices_range = current_index, [current_index]
    for o in self._order(order)[current_order:]:
      for p in xrange(o+1):
        for q in xrange(p+1): u += 1
      indices_range.append(u)

    # loop over remaining orders.
    for o, startu in zip(self._order(order)[current_order:], indices_range):
      # loop over results for each band.
      for band_x in all_xs.T:
        u = startu
        # create tensor from n vector.
        dummy = zeros([3]*o, dtype="float64")
        for p in xrange(o+1):
          for q in xrange(p+1):
            indices = [[i]*j for i, j in zip([0, 1, 2], [o-p, p-q, q]) if j != 0]
            indices = set(permutations(chain(*indices)))
            for index in indices: dummy[index] = band_x[u] if abs(band_x[u]) > tolerance else 0
            u += 1
        # add tensor to results for that order.
        result[current_order].append(dummy * self.eigenvalues.units * angstrom**o)
      current_order += 1

    # got it all.
    return result

  def inverse_mass_tensor(self, order=2, kmag=True, lstsq=None, tolerance=None, **kwargs):
    """ Returns inverse mass tensor, hopefully with right units. """
    from ..physics import electronic_mass, h_bar
    if 2 not in self._order(order):
      raise AttributeError("Effective mass are not part of the current expansion.")
    result = self.tensors(order, kmag, lstsq, tolerance, **kwargs)[self._order(order).index(2)]
    for i in xrange(len(result)): 
      result[i] = (result[i] / h_bar**2).rescale(1./electronic_mass)
    return result


    
def reciprocal( vasp, structure, outdir = None, comm = None,
                order = 2, nbpoints = None, stepsize = 1e-2, 
                center = None, eigtol = 1e-10, **kwargs ):
  """ Computes k-space taylor expansion of the eigenvalues up to given order.

      First runs a vasp calculation using the first input argument, regardless
      of whether a restart keyword argument is also passed. In practice,
      following the general LaDa philosophie of never overwritting previous
      calculations, this will not rerun a calculation if one exists in
      ``outdir``. 
      Second, a static non-self-consistent calculation is performed to compute
      the eigenvalues for all relevant kpoints.

      :Parameters:
        vasp : `vasp.Vasp` or `vasp.RelaxCellShape`
          Vasp or derived functional.
        structure : `crystal.Structure`
          The structure for wich to compute effective masses.
        outdir : str
          Root directory where to save results of calculation. Calculations
          will be stored in  "reciprocal" subdirectory of this input parameter.
        comm : `lada.mpi.Communicator` or None
          MPI communicator containing processes with which to perform
          calculation.
        order : int
          Highest order derivative to perform. Defaults to 1.
        nbpoints : int
          Number of points (in a single direction) with wich to compute
          taylor expansion.  Should be at least order + 1. Default to order +
          1. Note that in the case of even nbpoints, an extra kpoint is added
          so that the center does get computed.
        stepsize : float
          Distance between interpolation points. Default = 1e-2.
          Units of ``2|pi|/a``, with ``a=structure.scale``.
        center : 3-tuple
          Central k-point of the taylor expansion. This should be given in
          **reciprocal** units (eg coefficients to the reciprocal lattice
          vectors). Default is None and means |Gamma|.
        eigtol : float
          Energy convergence criteria (ediffg) for static calculation.
        kwargs 
          Extra parameters which are passed on first to escan (if escan
          object has an attribute with the same name), then to the linear least
          square fit method. Note that this will *not* change the external escan
          object.  This function is stateless. 

      :return: Extraction object from which derivatives can be obtained.

      .. |pi|     unicode:: U+003C0 .. GREEK SMALL LETTER PI
      .. |Gamma|  unicode:: U+00393 .. GREEK CAPITAL LETTER GAMMA
  """
  from copy import deepcopy
  from os import getcwd
  from os.path import join
  from numpy import array, dot, append
  from numpy.linalg import inv
  from . import Vasp
 
  # takes care of default parameters.
  if center is None: center = kwargs.pop("kpoint", [0,0,0])
  center = array(center, dtype="float64")
  if outdir is None: outdir = getcwd()
  maxorder = max(order) if hasattr(order, '__iter__') else order 
  if nbpoints == None: nbpoints = maxorder + 1
  if nbpoints < maxorder + 1:
    raise ValueError("Cannot compute taylor expansion of order {0} "\
                     "with only {1} points per direction.".format(maxorder, nbpoints))


  # first runs vasp.
  first = vasp.Extract(directory=outdir, comm=comm)
  if not first.success:
    first = vasp(structure=structure, outdir=outdir, comm=comm, **kwargs)
  if not first.success: return Extract(outdir=outdir, comm=comm)

  # prepare second run.
  functional = first.functional
  center = dot(inv(first.structure.cell).T, center)
  kpoints = [ (x, y, z) for x in xrange(nbpoints)\
                        for y in xrange(nbpoints)\
                        for z in xrange(nbpoints) ]
  functional.kpoints = (array(kpoints, dtype="float64") - (nbpoints-1)/2.) * stepsize + center
  if nbpoints % 2 == 0: # adds central point to grid.
    functional.kpoints = append(center[None,:], functional.kpoints, axis=0)

  # and exectute it.
  kwargs = deepcopy(kwargs)
  kwargs['restart']     = first
  kwargs['nonscf']      = True
  kwargs['relaxation']  = "static"
  kwargs['ediff']       = eigtol / len(structure.atoms)
  second = functional(first.structure, comm=comm, outdir=join(first.directory, "reciprocal"), **kwargs)
  return Extract(outcar=second.directory, comm=None, input=second)

reciprocal.Extract = Extract
""" Extractor class for the reciprocal method. """

  
