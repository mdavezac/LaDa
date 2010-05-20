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
    rspace2 = np.tensordot(v, extract.raw_wfns,[extract.inverse_indices,:,:].conjugate() ((1),(0)))
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


  
