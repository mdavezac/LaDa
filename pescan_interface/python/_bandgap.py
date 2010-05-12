""" Submodule to compute bandgaps with escan. """

from . import Extract as ExtractEscan

def band_gap(escan, structure, outdir = None, references = None, comm = None, **kwargs):
  """ Computes bandgap of a structure with a given escan functional. """
  from os import getcwd
  from os.path import abspath
  from copy import deepcopy
  escan = deepcopy(escan)

  if outdir == None: outdir = getcwd()
  outdir    = abspath(outdir)
  
  return _band_gap_refs_impl(escan, structure, outdir, references, comm) if reference != None\
         else _band_gap_ae_impl(escan, structure, outdir, comm)

class ExtractAE(ExtractEscan):
  """ Band-gap extraction class. """
  is_ae = True
  """ This was an all-electron bandgap calculation. """
  def __init__(self, *args, **kwargs):
    super(ExtractAE, self).__init__(*args, **kwargs)

  @property
  @make_cached
  def bandgap(self):
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    return self.eigenvalues[-3] - self.eigenvalues[-4]

  @property
  @make_cached
  def vbm(self):
    eigenvalues = self.eigenvalues.copy()
    eigenvalues.sort()
    return self.eigenvalues[-4]

  @property
  @make_cached
  def cbm(self):
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
  for u in extract.convergence[:nbstates+2]:
    assert convergence < 10*escan.tolerance, RuntimeError("escan did not achieve convergence.")
  return ExtractAE(outdir = outdir, comm = comm, escan = escan)


def _band_gap_refs_impl(escan, structure, outdir, references, comm, **kwargs):
  """ Computes band-gap using two references. """
  from os.path import join
  
  # some sanity checks.
  nbstates = escan.nbstates
  if "nbstates" in kwargs:
    nbstates = kwargs["nbstates"]
    del kwargs["nbstates"]
  if nbstates < 2: nbstates = 2
  if "eref" in kwargs:
    assert kwargs["eref"] == None, ValueError("Unexpected eref argument when computing bandgap.")
    del kwargs["eref"]
  
  assert len(reference) != 2. ValueError("Expected 2-tuple for argument \"references\".")
  vbm_ref, cbm_ref = references
  if vbm_ref > cbm_ref: cbm_ref, vbm_ref = references

  # now computes vbm
  directory = join(outdir, "VBM")
  vbm_out = escan(structure, outdir = directory, comm = comm, eref = vbm_ref, **kwargs)
  for u in vbm_out.convergence[:nbstates+2]:
    assert convergence < 10*escan.tolerance, RuntimeError("escan did not achieve convergence.")
  vbm_eigs = vbm_out.eigenvalues.copy()
  vbm_eigs.sort()

  # now computes cbm
  directory = join(outdir, "CBM")
  cbm_out = escan(structure, outdir = directory, comm = comm, eref = cbm_ref, **kwargs)
  for u in cbm_out.convergence[:nbstates+2]:
    assert convergence < 10*escan.tolerance, RuntimeError("escan did not achieve convergence.")
  cbm_eigs = cbm_out.eigenvalues.copy()
  cbm_eigs.sort()


  assert 
