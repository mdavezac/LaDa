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

def scale(structure):
  """ Approximate scale of a given structure. """
  n = len(structure.atoms)
  nSi = len([a.type for a in structure.atoms if a.type == 'Si'])
  return float(nSi) / float(n) * structure.scale + float(n - nSi) / float(n) * 5.69

def direction001(n):
  """" Cell and direction for 001. """
  a, b = (n/2, 0) if n % 2 == 0 else (n/2+0.5, 0.5)
  return [ [a, 0, 0], [b, 0.5, 0.5], [0, -0.5, 0.5] ], [1, 0, 0]
def direction011(n):
  """" Cell and direction for 011. """
  return [[0, 0.5, 0.5],[n*0.5,0, 0.5], [n*0.5,0.5,0]], [0, 1, 1]
def direction111(n):
  """" Cell and direction for 111. """
  return [[0.5*n, -0.5, 0.5],[0.5*n,0, -0.5], [0.5,0.5,0]], [1, 1, 1]

def create_sl(path, direction, nmin, nmax, nstep, x=0.5, input='input.py'):
  """ Creates dictionary with input structures. """
  from IPython.ipapi import get as get_ipy
  from numpy import array
  from numpy.linalg import norm, inv
  from lada.jobs import JobDict
  from lada.escan import read_input, exec_input, ReducedKDensity
  from lada.crystal.binary import zinc_blende
  from lada.crystal import nb_valence_states, layer_iterator

  ip = get_ipy()

  input = read_input(input)
  kescan = exec_input(repr(input.escan).replace('Escan', 'KEscan')).functional

  lattice = zinc_blende()
  lattice.sites[0].type = 'Si', 'Ge'
  lattice.sites[1].type = 'Si', 'Ge'
  lattice.scale = 5.45
  lattice.find_space_group()

  density = 10e0 * max([1e0/norm(u) for u in inv(lattice.cell)])

  rootjobs = ip.user_ns.get('current_jobdict', JobDict())
  for n0 in range(nmin, nmax, nstep):
    # create structure
    n1 = int(n0 / x - n0 + 0.5)
    cell, dir = direction(n0+n1)
    structure = lattice.to_structure(cell)
    N0 = len(structure.atoms) * x
    for layer in layer_iterator(structure, dir):
      for atom in layer:
        atom.type = 'Si' if N0 > 0 else 'Ge'
        N0 -= 1
    # name and scale.
    structure.name = "{0[0]}{0[1]}{0[2]}/x_{1:0<4.3}/n_{2}".format(dir, x, n0)
    structure.scale = scale(structure)
    # creates jobdictionary.
    jobdict = rootjobs / structure.name
    jobdict.jobparams['structure'] = structure.copy()
    jobdict.functional = kescan.copy()
    jobdict.functional.kpoints = ReducedKDensity(density, (0.5, 0.5, 0.5))
    jobdict.functional.reference = None
    jobdict.functional.fft_mesh = fftmesh(structure.cell)
    jobdict.functional.nbstates = int(nb_valence_states(structure) * 1.5+0.5)

  if 'current_jobdict' not in ip.user_ns: ip.user_ns["current_jobdict"] = rootjobs
  ip.magic("savejobs " + path)
  return
