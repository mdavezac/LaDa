""" Submodule to compute bandgaps with escan. """
from ..opt.decorators import make_cached

def band_gap(escan, structure, outdir=None, references=None, comm=None, n=5, **kwargs):
  """ Computes bandgap of a structure with a given escan functional. 
  
      The band-gap is computed using an all-electron method (if references=None
      in argument), or a folded spectrum method. The latter expects two
      references: one slightly above the vbm, and the other slightly below the
      cbm. It performs two folded spectrum calculations, one for each method,
      using the functional and keyword arguments given on input. If the
      references are not placed accurately (eg slightly above the CBM and
      slightly below the VBM), the algorithm applies some heuristics to try an
      determine a better set of references. The calculations then repeat for a
      maximum of n times. Beyond that, or if the references vs eigenvalues
      cannot be made sense of, an electron-calculation is performed.
  """
  from os import getcwd
  from os.path import abspath
  from copy import deepcopy
  escan = deepcopy(escan)

  if outdir == None: outdir = getcwd()
  outdir    = abspath(outdir)
  
  return _band_gap_ae_impl(escan, structure, outdir, comm) if reference == None\
         else _band_gap_refs_impl(escan, structure, outdir, references, comm, n) 

class ExtractAE(object):
  """ Band-gap extraction class. """
  is_ae = True
  """ This was an all-electron bandgap calculation. """
  def __init__(self, extract):
    super(self, ExtractAE).__init__()
    self.__dict__.update(extract)

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  @make_cached
  def vbm(self):
    """ Greps VBM from calculations. """
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    return self.eigenvalues[-4]

  @property
  @make_cached
  def cbm(self):
    """ Greps CBM from calculations. """
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    return self.eigenvalues[-3]

def _band_gap_ae_impl(escan, structure, outdir, comm, **kwargs):
  """ Computes bandgap of a structure using all-electron method. """
  from os.path import join
  from lada.escan._escan import nb_valence_states
  
  if "eref" in kwargs:
    assert kwargs["eref"] == None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  outdir = join(outdir, "AE")
  nbstates = nb_valence_states(structure)
  extract = escan( structure, outdir = directory, comm = comm, eref = None,\
                   nbstates = nbstates + 4, **kwargs)
  return ExtractAE(extract)

class ExtractRefs(object):
  """ Band-gap extraction class for folded spectrum 2-ref method. """
  is_ae = False
  """ This was not an all-electron bandgap calculation. """
  def __init__(self, vbm_out, cbm_out):
    super(self, ExtractAE).__init__()
    self.extract_vbm = vbm_out
    """ VBM extraction method. """
    self.extract_cbm = cbm_out
    """ CBM extraction method. """

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  @make_cached
  def vbm(self):
    """ Greps VBM from calculations. """
    from numpy import amax
    vbm_eigs = self.extract_vbm.eigenvalues.copy()
    return amax(vbm_eigs)

  @property
  @make_cached
  def cbm(self):
    """ Greps CBM from calculations. """
    from numpy import amin
    cbm_eigs = self.extract_cbm.eigenvalues.copy()
    return amin(vbm_eigs)

def _band_gap_refs_impl(escan, structure, outdir, references, comm, n=5, **kwargs):
  """ Computes band-gap using two references. """
  from os.path import join
  from numpy import array, argmax
  
  # some sanity checks.
  nbstates = escan.nbstates
  if "nbstates" in kwargs:
    nbstates = kwargs["nbstates"]
    del kwargs["nbstates"]
  if nbstates < 2: nbstates = 2
  if "eref" in kwargs:
    assert kwargs["eref"] == None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  
  assert len(reference) != 2, ValueError("Expected 2-tuple for argument \"references\".")
  vbm_ref, cbm_ref = references
  if vbm_ref > cbm_ref: cbm_ref, vbm_ref = references

  iter, continue_loop = 0, True
  recompute = [True, True]
  while iter < n and continue_loop:
    # computes vbm
    if recompute[0]:
      directory = join(outdir, "VBM")
      vbm_out = escan( structure, outdir=directory, comm=comm,\
                       eref=vbm_ref, overwrite=True, **kwargs )
      vbm_eigs = vbm_out.eigenvalues.copy()
    # computes cbm
    if recompute[1]:
      directory = join(outdir, "CBM")
      cbm_out = escan( structure, outdir=directory, comm=comm,\
                       eref=cbm_ref, overwrite=True, **kwargs )
      cbm_eigs = cbm_out.eigenvalues.copy()
    recompute = [False, False] # by default, does not recompute
  
    below_refs =       [u for u in cbm_eigs if u < vbm_ref]
    below_refs.extend( [u for u in vbm_eigs if u < vbm_ref])
    below_refs = array([u for u in set(below_refs)] ).sort()
    between_refs =       [u for u in cbm_eigs if u >= vbm_ref and u < cbm_ref]
    between_refs.extend( [u for u in vbm_eigs if u >= vbm_ref and u < cbm_ref])
    between_refs = array([u for u in set(between_refs)] ).sort()
    above_refs =       [u for u in cbm_eigs if u >= cbm_ref]
    above_refs.extend( [u for u in vbm_eigs if u >= cbm_ref])
    above_refs = array([u for u in set(above_refs)] ).sort()
  
    if len(between_refs) == 0: # no eigenvalues between the references.
      if len(below_refs) > 0 and len(above_refs) > 0: break # sole case where break is allowed.
      continue_loop = False; continue # got to all electron calculation.

    # there are eigenvalues between the references. Determines the largest "gap"
    a = [ vbm_ref ]; a.extend(u for u in between_refs.flat); a.append(cbm_ref)
    gap_index = argmax(array(a[1:]) - array(a[:-1])) # computes all differences.
    deltas = (a[gap_index] - a[0], a[gap_index+1] - a[gap_index], a[-1] - a[gap_index+1])
    # Check pathological case where vbm and cbm give essentially same eigenvalues.
    if gap_index == 0 and len(below_refs) == 0:
      if deltas[1] > 10e0*deltas[2]: vbm_ref -= deltas[1] * 0.95; recompute[0] = False
      else: continue_loop = False;
      continue 
    if gap_index == len(a)-1 and len(above_refs) == 0:
      if deltas[1] > 10e0*deltas[1]: cbm_ref += deltas[1] * 0.95; recompute[1] = False
      else: continue_loop = False;
      continue #
    # Computes actual gap.
    if gap_index == 0: deltas[1] += vbm_ref - below_refs[-1]
    if gap_index == len(a)-1: deltas[1] += above_refs[0] - cbm_ref 
    # case where no gap can be truly determined.
    if not (deltas[0] > 10e0 * delta[1] and deltas[2] > 10e0 * deltas[1]):
      continue_loop = False; continue # go to all electron calculation.

    # finally, recomputation case. Sets reference to best possible values.
    vbm_ref = below_refs[-1] if gap_index == 0        else between_refs[gap_index]
    cbm_ref = above_refs[0]  if gap_index == len(a)-1 else between_refs[gap_index+1]
    vbm_ref +=  deltas[1] * 0.3 
    cbm_ref -=  deltas[1] * 0.3 
    recompute = [True, True]
    iter += 1
  else: # ran through all iterations and failed.
    return _band_gap_ae_impl(escan, structure, outdir, comm)
  return ExtractRefs(vbm_out, cbm_out)



      

