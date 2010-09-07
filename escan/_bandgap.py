""" Submodule to compute bandgaps with escan. """
__docformat__  = 'restructuredtext en'
from ..opt.decorators import make_cached
from _extract import Extract as _ExtractE

def extract(outdir=".", comm = None):
  """ Gets extraction object from directory structure. 
  
      Checks first for all-electron calculation in outdir/AE. If it does not
      exist or the calculation is unsucessful, then checks for VBM, CBM
      directories, and checks those calculations. If unsucessfull, returns a
      fake extraction object with success set to False. Otherwise, returns
      successfull extraction object.
  """
  from os.path import exists, join
  from ..vff import Extract as VffExtract
  from boost.mpi import broadcast

  is_root = True if comm == None else comm.rank == 0

  paths = join(outdir, "AE"), join(outdir, "VBM"), join(outdir, "CBM")
  if is_root: exists_paths = [exists(p) for p in paths]
  if comm != None:
    exists_paths = broadcast(comm, exists_paths if is_root else None, root=0)
  if exists_paths[0]: 
    result = ExtractAE( _ExtractE(paths[0], comm = comm) )
    if result.success: return result
  elif exists_paths[1] and exists_paths[2]:
    result = ExtractRefs( _ExtractE(paths[1], comm = comm),\
                          _ExtractE(paths[2], comm = comm),
                          _ExtractE(outdir, comm = comm) )
    if result.success: return result
  class NoBandGap(object):
    def __init__(self): self.success = False
  return NoBandGap()

def band_gap(escan, structure, outdir=None, references=None, n=5, overwrite = False, **kwargs):
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
  from os.path import abspath, exists, join
  from copy import deepcopy
  from boost.mpi import world, broadcast

  assert "do_escan" not in kwargs,\
         ValueError("\"do_escan\" is not an admissible argument of bandgap.")

  escan = deepcopy(escan)
         
  if outdir == None: outdir = getcwd()
  outdir    = abspath(outdir)
  overlap_factor = kwargs.pop("overlap_factor", 10e0)

  comm = kwargs.pop("comm", world)
  if not overwrite:  # check for previous results.
    result = extract(outdir, comm)
    if result.success: return result

  
  kwargs["overwrite"] = overwrite
  kwargs["comm"] = comm
  return _band_gap_ae_impl(escan, structure, outdir, **kwargs) if references == None\
         else _band_gap_refs_impl(escan, structure, outdir, references, n, \
                                  overlap_factor=overlap_factor, **kwargs) 

class ExtractAE(_ExtractE):
  """ Band-gap extraction class. """
  is_ae = True
  """ This was an all-electron bandgap calculation. """
  def __init__(self, extract):
    super(ExtractAE, self).__init__(directory=extract.directory, comm=extract.comm)
    self.OUTCAR = extract.OUTCAR

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  def vbm(self):
    """ Greps VBM from calculations. """
    return self._vbm_cbm[0]

  @property
  def cbm(self):
    """ Greps CBM from calculations. """
    return self._vbm_cbm[1]

  @property
  @make_cached
  def _vbm_cbm(self):
    """ Gets vbm and cbm. """
    from quantities import eV
    from lada.escan._escan import nb_valence_states
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    n = nb_valence_states(self.structure)
    a, b = eigenvalues[n-1], eigenvalues[n] 
    if hasattr(a, "units"): a.units = eV
    else: a = a * eV
    if hasattr(b, "units"): b.units = eV
    else: b = b * eV
    return a, b


  def oscillator_strength(self, degeneracy=1e-3, attenuate=False):
    """ Computes oscillator strength between vbm and cbm. """
    from numpy import dot
    from numpy.linalg import det
    from ..physics import a0, electronic_mass, h_bar
    from . import dipole_matrix_elements
    result, nstates = None, 0
    units = 2e0/3e0 * h_bar**2 / electronic_mass
    for wfnA in self.gwfns:
      if abs(wfnA.eigenvalue - self.cbm) > degeneracy: continue
      for wfnB in self.gwfns:
        if abs(wfnB.eigenvalue - self.vbm) > degeneracy: continue
        nstates += 1
        dme = wfnA.braket(self.gvectors, wfnB, attenuate=attenuate) 
        if result == None: 
          result = dot(dme, dme.conjugate()).real / (wfnA.eigenvalue - wfnB.eigenvalue) \
                   * dme.units * dme.units
        else: 
          result += dot(dme, dme.conjugate()).real / (wfnA.eigenvalue - wfnB.eigenvalue) \
                    * dme.units * dme.units
    return (result * units).simplified, nstates
  
def _band_gap_ae_impl(escan, structure, outdir, **kwargs):
  """ Computes bandgap of a structure using all-electron method. """
  from os.path import join
  from lada.escan._escan import nb_valence_states
  
  if "eref" in kwargs:
    assert kwargs["eref"] == None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  outdir = join(outdir, "AE")
  nbstates = kwargs.pop("nbstates", escan.nbstates)
  if nbstates == None: nbstates = 4
  if nbstates == 0: nbstates = 4
  nbstates = nbstates  + nb_valence_states(structure)
  extract = escan( structure, outdir = outdir, eref = None,\
                   nbstates = nbstates, **kwargs)
  return ExtractAE(extract)

class ExtractRefs(object):
  """ Band-gap extraction class for folded spectrum 2-ref method. """
  is_ae = False
  """ This was not an all-electron bandgap calculation. """
  def __init__(self, vbm_out, cbm_out, vff_out):
    super(ExtractRefs, self).__init__()
    self.extract_vbm = vbm_out
    """ VBM extraction method. """
    self.extract_cbm = cbm_out
    """ CBM extraction method. """
    self.extract_vff = vff_out
    """ VFF extraction method. """

  @property
  def bandgap(self):
    """ Greps band-gap from calculation. """
    return self.cbm - self.vbm

  @property
  @make_cached
  def _raw(self):
    from numpy import array
    from quantities import eV
    vbm_eigs = self.extract_vbm.eigenvalues.copy()
    cbm_eigs = self.extract_cbm.eigenvalues.copy()

    vbm_ref = self.extract_vbm.escan.eref
    cbm_ref = self.extract_cbm.escan.eref
    above_refs =       [u for u in cbm_eigs if u >= cbm_ref]
    above_refs.extend( [u for u in vbm_eigs if u >= cbm_ref])
    above_refs = array([u for u in set(above_refs)] );
    above_refs.sort()
    below_refs =       [u for u in cbm_eigs if u <= vbm_ref]
    below_refs.extend( [u for u in vbm_eigs if u <= vbm_ref])
    below_refs = array([u for u in set(below_refs)] );
    below_refs.sort()

    a, b = below_refs[-1], above_refs[0]
    if hasattr(a, "units"): a.units = eV
    else: a = a * eV
    if hasattr(b, "units"): b.units = eV
    else: b = b * eV
    return a, b

  @property
  def vbm(self):
    """ Greps VBM from calculations. """
    from numpy import amax
    vbm_eigs = self.extract_vbm.eigenvalues.copy()
    return self._raw[0]

  @property
  def cbm(self):
    """ Greps CBM from calculations. """
    from numpy import amin
    cbm_eigs = self.extract_cbm.eigenvalues.copy()
    return self._raw[1]

  def oscillator_strength(self, degeneracy=1e-3, attenuate=False):
    """ Computes oscillator strength between vbm and cbm. """
    from numpy import all, abs, dot
    from numpy.linalg import det
    from ..physics import a0, electronic_mass, h_bar
    from . import dipole_matrix_elements

    assert self.extract_vbm.comm == self.extract_cbm.comm
    assert self.extract_vbm.gvectors.shape == self.extract_cbm.gvectors.shape
    assert all( abs(self.extract_vbm.gvectors - self.extract_cbm.gvectors) < 1e-12 )
    units = 2e0/3e0 * h_bar**2 / electronic_mass
    result, nstates = None, 0
    for wfnA in self.extract_cbm.gwfns:
      if abs(wfnA.eigenvalue - self.cbm) > degeneracy: continue
      for wfnB in self.extract_vbm.gwfns:
        if abs(wfnB.eigenvalue - self.vbm) > degeneracy: continue
        nstates += 1
        dme = wfnA.braket(self.extract_vbm.gvectors, wfnB, attenuate=attenuate)
        dme = dot(dme, dme.conjugate()).real * dme.units * dme.units
        if result == None: result = dme / (wfnA.eigenvalue - wfnB.eigenvalue) 
        else: result += dme / (wfnA.eigenvalue - wfnB.eigenvalue) 
    return (units * result).simplified, nstates

  @property
  def success(self):
    """ True if all calculations were successfull. """
    return self.extract_vff.success and self.extract_vbm.success and self.extract_cbm.success

  def __getattr__(self, name): 
    """ Sends to vff output. """
    if hasattr(self.extract_vff, name): return getattr(self.extract_vff, name)
    raise AttributeError("Unknown attribute %s." % (name))

  def __dir__(self):
    """ Returns list of attributes. 

        Since __getattr__ is modified to include vff extraction data, __dir__
        should be modified as well. This makes for better ipython
        (tab-completion) integration. 
    """
    result = [u for u in dir(self.__class__) if u[0] != '_'] 
    result.extend(['extact_vff', 'extract_vbm', 'extract_cbm', 'structure', 'escan', 'vff'])
    return list(set(result))

def _band_gap_refs_impl( escan, structure, outdir, references, n=5,\
                         overlap_factor=10e0, **kwargs):
  """ Computes band-gap using two references. """
  from os.path import join, exists
  from shutil import copyfile
  from numpy import array, argmax
  
  # check/correct input arguments
  if "overwrite" in kwargs: del kwargs["overwrite"]
  nbstates = kwargs.pop("nbstates", escan.nbstates)
  if nbstates < 2: nbstates = 2
  if "eref" in kwargs:
    assert kwargs["eref"] == None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  vffrun = kwargs.pop("vffrun", escan.vffrun)
  genpotrun = kwargs.pop("genpotrun", escan.genpotrun)
  
  assert len(references) == 2, ValueError("Expected 2-tuple for argument \"references\".")
  vbm_ref, cbm_ref = references
  if vbm_ref > cbm_ref: cbm_ref, vbm_ref = references

  # first computes vff and genpot unless given.
  if genpotrun == None or vffrun == None: 
    vffout = escan( structure, outdir=outdir, do_escan=False, genpotrun=genpotrun,\
                    vffrun=vffrun, **kwargs )
    if genpotrun == None: genpotrun = vffout
    if vffrun == None: vffrun = vffout
  else: vffout = vffrun

  iter, continue_loop = 0, True
  recompute = [True, True]
  while iter < n and continue_loop:
    # computes vbm
    if recompute[0]:
      vbm_out = escan\
                (
                  structure, outdir=join(outdir,"VBM"), \
                  eref=vbm_ref, overwrite=True, vffrun=vffrun,\
                  genpotrun=genpotrun, nbstates=nbstates, **kwargs 
                )
      vbm_eigs = vbm_out.eigenvalues.copy()
    # computes cbm
    if recompute[1]:
      cbm_out = escan\
                (
                  structure, outdir=join(outdir, "CBM"), \
                  eref=cbm_ref, overwrite=True, vffrun=vffrun,
                  genpotrun=genpotrun, nbstates=nbstates, **kwargs
                )
      cbm_eigs = cbm_out.eigenvalues.copy()
    recompute = [False, False] # by default, does not recompute
  
    below_refs =       [u for u in cbm_eigs if u < vbm_ref]
    below_refs.extend( [u for u in vbm_eigs if u < vbm_ref])
    below_refs = array([u for u in set(below_refs)] )
    below_refs.sort()
    between_refs =       [u for u in cbm_eigs if u >= vbm_ref and u < cbm_ref]
    between_refs.extend( [u for u in vbm_eigs if u >= vbm_ref and u < cbm_ref])
    between_refs = array([u for u in set(between_refs)] );
    between_refs.sort()
    above_refs =       [u for u in cbm_eigs if u >= cbm_ref]
    above_refs.extend( [u for u in vbm_eigs if u >= cbm_ref])
    above_refs = array([u for u in set(above_refs)] )
    above_refs.sort()
  
    if between_refs.size == 0: # no eigenvalues between the references.
      if below_refs.size > 0 and above_refs.size > 0: break # sole case where break is allowed.
      continue_loop = False; continue # got to all electron calculation.

    # there are eigenvalues between the references. Determines the largest "gap"
    a = [ vbm_ref ]; a.extend(u for u in between_refs.flat); a.append(cbm_ref)
    gap_index = argmax(array(a[1:]) - array(a[:-1])) # computes all differences.
    deltas = [a[gap_index] - a[0], a[gap_index+1] - a[gap_index], a[-1] - a[gap_index+1]]
    # Check pathological case where vbm and cbm give essentially same eigenvalues.
    if gap_index == 0 and below_refs.size == 0:
      if deltas[1] > overlap_factor * deltas[2]: vbm_ref -= deltas[1] * 0.95; recompute[0] = False
      else: continue_loop = False;
      continue 
    if gap_index == len(a)-1 and above_refs.size == 0:
      if deltas[1] > overlap_factor * deltas[1]: cbm_ref += deltas[1] * 0.95; recompute[1] = False
      else: continue_loop = False;
      continue #
    # Computes actual gap.
    if gap_index == 0: deltas[1] += vbm_ref - below_refs[-1]
    if gap_index == len(a)-1: deltas[1] += above_refs[0] - cbm_ref 
    # case where no gap can be truly determined.
    if not (deltas[1] > overlap_factor * deltas[0] and deltas[1] > overlap_factor * deltas[2]):
      continue_loop = False; continue # go to all electron calculation.

    # finally, recomputation case. Sets reference to best possible values.
    vbm_ref = below_refs[-1] if gap_index == 0        else a[gap_index]
    cbm_ref = above_refs[0]  if gap_index == len(a)-1 else a[gap_index+1]
    vbm_ref +=  deltas[1] * 0.3 
    cbm_ref -=  deltas[1] * 0.3 
    recompute = [True, True]
    iter += 1
  else: # ran through all iterations and failed.
    return _band_gap_ae_impl(escan, structure, outdir, **kwargs)
  return ExtractRefs(vbm_out, cbm_out, vffout)






      

