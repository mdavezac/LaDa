from lada.enumeration import Enum
from lada.crystal.binary import zinc_blende
from lada.jobs import JobDict
from lada.escan import read_input, exec_input

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

input = read_input('input.py')
kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

enum = Enum(zinc_blende())
enum.sites[0].type = 'Si', 'Ge'
enum.sites[1].type = 'Si', 'Ge'
enum.find_space_group()

jobs = JobDict()

i = 0
for n in range(1, 3):
  for structure in enum.structures(n): pass
      structure.name = str(i)
      nSi = len([a.type for a in structure.atoms if a.type == 'Si'])
      structure.scale = float(nSi) / float(n) * 5.45 + float(n - nSi) / float(n) * 5.65

      jobdict = jobs / structure.name
      jobdict.jobparams['structure'] = structure.copy()
      jobdict.jobparams['fftmesh'] = fftmesh(structure.cell)
      jobdict.functional = kescan
      print structure

      i += 1
  print i
