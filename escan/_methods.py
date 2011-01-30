""" Subpackage containing dipole matrix-element implementations. """

def dipole_matrix_elements(*args):
  """ Computes dipole matrix elements (g-space).

      If only one extraction/output object is given, the dipole moments are
      computed between all states in that particular calculation.  If two
      extraction/output objects are given, then dipole matrix elements between
      the wavefunctions in those calculations are given (eg inter, not intra).
      All calculations are performed in g-space, within the approximation that
      [r, H] ~ ip.
  """
  from numpy import zeros, all, transpose

  if len(args) == 1:  outA, outB = args[0], args[0]
  elif len(args) == 2:
    outA, outB = args[0], args[1]
    assert all(abs(outA.gvectors-outB.gvectors) < 1e-8),\
           RuntimeError("Calculations were not performed within the same conditions.")
  else: raise RuntimeError("Unknown arguments to dipole_matrix_elements.")

  result = zeros((outA.escan.nbstates, outB.escan.nbstates, 3), dtype="float64")
  
  for i, awfn in enumerate(outA.gwfns):
    for j, bwfn in enumerate(outB.gwfns):
      result[i,j,:] = awfn.braket(transpose(outA.gvectors), bwfn, attenuate=True)

  return -1e0j * result


def LDOS(extract, positions):
  """ Local density of states """
  import numpy as np
  from boost.mpi import broadcast, reduce

  result = None
  # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
  v = np.exp(-1j * np.tensordot(positions, extract.gvectors, ((1),(1))))
  # computes fourrier transform for all wavefunctions simultaneously.
  rspace = np.tensordot(v, extract.raw_wfns, ((1),(0)))
  # reduce across processes
  rspace = reduce(comm, result, lambda x,y: x+y, 0)
  if extract.is_krammer:
    rspace2 = rspace 
    rspace = zeros( (rspace.shape[0], rspace.shape[1]*2, rspace.shape[2]), dtype="complex64")
    rspace[:,::2,:] = rspace2
    rspace2 = np.tensordot(v, extract.raw_wfns[extract.inverse_indices,:,:].conjugate(), ((1),(0)))
    rspace2 = reduce(comm, result, lambda x,y: x+y, 0)
    rspace[:,1::2,:] = rspace2
  rspace = np.multiply(rspace, np.conjugate(rspace))
  if not extract.is_spinor: rspace = rspace[:,:,0]
  else: rspace = rspace[:,:,0] + rspace[:,:,1]

 
  class ldoses(object):
    """ Local density of states for a given set of positions within a given structure. """
    def __init__(self, eigs, rs):
      self.eigs, self.rs = -np.multiply(eigs, eigs), rs.copy()
    def __call__(self, e, sigma=0.1):
      return np.dot(self.rs, 1e0/sqrt(np.pi)/sigma * np.exp(-self.eigs/sigma/sigma))

  return ldoses(extract.eigenvalues, rspace)

def majority_representation(extractor, multicell, tolerance=1e-12): 
  """ Majority Represention for all unmapped kpoint.
  
      :Parameters: 
        extractor : `KExtractor`
          Result of a `KEscan` calculation. The ``kpoints`` attribute from
          ``extractor.functional`` should be a class of *reduced* kpoints, such
          as `ReducedKGrid`, `ReducedKDensity`, or `ReducedBPoints`.
        multicell : numpy 3x3
          A matrix (M) linking the unit-cell (U) to the super cell (S), ``S = U
          * M``.  From this, the unit-cell is determined, whatever deformation
          may have happened, as engendered by VFF.
       
       :return: An  list of lists of 2-tuples. The first dimension corresponds to unreduced kpoints.
         The 2-tuples consiste  of the eigenvalue for that kpoint and band, and
         of its majority representation. The inner list duplicate entries with
         respect to eigenvalues. It is meant to be used with matplotlib's
         hexbin.
  """
  from operator import itemgetter
  from numpy import dot, pi, array
  from numpy.linalg import inv
  from quantities import angstrom
  from ..crystal import is_on_lattice, to_voronoi
  from ..physics import a0

  assert extractor.success,\
         ValueError("extractor is not the return of a *successfull* KEscan calculation.")
  assert len(extractor) > 0, ValueError("extractor does not contain escan runs.")
  istr = extractor.input_structure
  ostr = extractor.structure
  unitcell = dot(ostr.cell, inv(multicell)) * ostr.scale * angstrom
  unitkcell = inv(unitcell.T) * 2e0 * pi / angstrom
  invcell = inv(ostr.cell.T)
  
  # use the functional's kpoint object so we have access to all its goodies.
  kpoints = extractor.functional.kpoints

  # list of results.
  results = []

  unreduced = array([k[1] for k in kpoints.unreduced(istr, ostr)])
  mapping = [i for i in kpoints.mapping(istr, ostr)]
  scale = float(1e0/(ostr.scale  * angstrom).rescale(a0))
  N = float(len(ostr.atoms))
  for k0, extract in zip(unreduced, extractor):
#   # actual calculation for that kpoint.
#   extract = extractor[i]
    # reciprocal vectors of supercell.  
    Gs = extract.gvectors        
    # Compute kpoint folded into voronoi cell.
    K = to_voronoi(k0, invcell)
    # links reduced to unreduced: folding vector.
    Gfold = (K - k0) * scale
    # mapping for reciprocal vector of supercell which are reciprocal vectors
    # of unitcell *with* folding.
    onlatop = is_on_lattice(Gs, unitkcell.rescale(Gs.units), Gfold)
    # loop over bands and sum 
    current = [[w.eigenvalue, w.expectation_value(onlatop).real] for w in extract.gwfns]
    current = sorted(current, key=itemgetter(0))
    intermediate = [current[0]]
    for eig, val in current[1:]: 
      if abs(eig-intermediate[-1][0]) < tolerance: intermediate[-1][1] += val
      else: intermediate.append([eig, val])

    results.append([ (e, v) for e, v in intermediate if abs(v) * N > tolerance])

  return results
      


  
