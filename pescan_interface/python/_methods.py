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
  from numpy import zeros, all

  if len(args) == 1:  outA, outB = args[0], args[0]
  elif len(args) == 2:
    outA, outB = args[0], args[1]
    assert all(abs(outA.gvectors-outB.gvectors) < 1e-8),\
           RuntimeError("Calculations were not performed within the same conditions.")
  else: raise RuntimeError("Unknown arguments to dipole_matrix_elements.")

  result = zeros((outA.escan.nbstates, outB.escan.nbstates, 3), dtype="float64")
  
  for i, awfn in enumerate(outA.gwfns):
    for j, bwfn in enumerate(outB.gwfns):
      result[i,j,:] = awfn.braket(outA.gvectors, bwfn, attenuate=True)

  return -1e0j * result

