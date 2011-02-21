def fftmesh(cell, cutoff=20e0):
  """ FFT real-space mesh used by escan. """
  from numpy.linalg import norm
  result = int(norm(cell[:,0]) * cutoff + 5e-1), \
           int(norm(cell[:,1]) * cutoff + 5e-1), \
           int(norm(cell[:,2]) * cutoff + 5e-1)
  # This line makes sure mesh is even. Not sure why...
  result = result[0] + result[0] % 2, \
           result[1] + result[1] % 2, \
           result[2] + result[2] % 2
  return result

def create_start(path, input='input.py', nall = 3, nrand = 5, nmax=100):
  """ Creates dictionary with input structures. """
  from random import shuffle
  from math import sqrt
  from itertools import chain
  from IPython.ipapi import get as get_ipy
  from lada.enumeration import Enum
  from lada.crystal.binary import zinc_blende
  from lada.jobs import JobDict
  from lada.escan import read_input, exec_input, ReducedKDensity

  input = read_input(input)
  kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

  enum = Enum(zinc_blende())
  enum.sites[0].type = 'Si', 'Ge'
  enum.sites[1].type = 'Si', 'Ge'
  enum.find_space_group()

  strs = [u for  n in range(nall, nrand) for u in enum.xn(n)]
  shuffle(strs)
  strs = [enum.as_structure(*u) for u in strs[:nmax]]
  alls = [structure for n in range(nall) for structure in enum.structures(n)]

  jobs = JobDict()
  for i, structure in enumerate(chain(alls, strs)):
    structure.name = str(i)
    nSi = len([a.type for a in structure.atoms if a.type == 'Si'])
    structure.scale = float(nSi) / float(n) * 5.45 + float(n - nSi) / float(n) * 5.65

    jobdict = jobs / structure.name
    jobdict.jobparams['structure'] = structure.copy()
    jobdict.jobparams['fftmesh'] = fftmesh(structure.cell)
    jobdict.functional = kescan
    jobdict.functional.kpoints = ReducedKDensity(1e0/20e0/sqrt(2e0), 0.5)
    jobdict.functional.reference = None

  ip = get_ipy()
  ip.user_ns["current_jobdict"] = jobdict.root
  ip.magic("savejobs " + path)
  return
