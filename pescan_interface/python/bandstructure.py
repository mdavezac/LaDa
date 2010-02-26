#! python
""" Bandstructure plotting tools """

def _line(start, end, density):
  """ Generator for creating points between two kpoints. """
  from numpy.linalg import norm

  distance = norm(end - start) 
  nbkpt = int(max(1, density * distance - 1))
  stepsize = 1e0/float(nbkpt)
  kpoints = [ float(i) * stepsize for i in range(1, nbkpt+1) ]
  for k in kpoints: yield start + k * (end-start) 

def _lines(endpoints, density):
  """ Generator for creating segments. """
  from numpy.linalg import norm
  
  assert len(endpoints) > 0, ValueError
  assert len(endpoints[0]) == 2, ValueError

  pos = 0
  yield pos, endpoints[0][0]
  for start, end in endpoints:
    for kpoint in _line(start, end, density):
      yield pos+norm(kpoint-start), kpoint
    pos += norm(end-start) 

  

def band_structure(structure, escan, vff, kpoints, density, **kwargs ):
  """ Returns eigenvalues for plotting bandstructure. """
  from os.path import join
  from copy import deepcopy
  from math import pow, pi, factorial
  from numpy import zeros, array, dot
  from numpy.linalg import norm, lstsq as np_lstsq
  from ._escan import method, potential
  from ..physics import a0, Hartree

  # saves mpicomm if necessary.
  mpicomm = escan.mpicomm
  if "mpicomm" in kwargs:
    mpicomm = kwargs["mpicomm"] 
    del kwargs["mpicomm"]
  # now copies escan
  escan = deepcopy(escan)
  # resets mpicomm to whatever it should be
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
  nbstates = escan.nbstates

  results = []
  for i, (x, escan.kpoint) in enumerate(_lines(kpoints, density)):
    # checks for double/krammer mad degeneracy touble
    double_trouble = escan.potential != potential.spinorbit or norm(escan.kpoint) < 1e-12
    # in which case only half the eigenvalues are computed.
    if double_trouble: escan.nbstates = nbstates / 2
    # sets directory.
    escan.directory = join(join(directory, "band_structure"), "%i-%s" % (i, escan.kpoint))
    # actually computes stuff.
    eigenvalues = escan(vff, structure)
    eigenvalues.sort()
    if double_trouble: # in case escan tries to screw us up again.
      eigenvalues = array([u for j in range(2) for u in eigenvalues])
      escan.nbstates = nbstates
    results.append( (x, escan.kpoint, eigenvalues) )
  return results
