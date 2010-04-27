
def create_zb_lattice():
  from numpy import array as np_array
  from lada.crystal import Lattice, Site

  lattice = Lattice()
  lattice.cell = np_array( [ [  0, 0.5, 0.5],
                             [0.5,   0, 0.5],
                             [0.5, 0.5,   0] ] )
  lattice.sites.append( Site( np_array([0.0,0,0]), ["Si", "Ge"] ) )
  lattice.sites.append( Site( np_array([0.25,0.25,0.25]), ["Si", "Ge"] ) )
  lattice.scale = 5.45
  lattice.find_space_group()
  lattice.set_as_crystal_lattice() 
  return lattice

def gtor_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  import numpy as np
  from boost.mpi import broadcast, reduce

  result = None
  for node in range(comm.size):
    # sends rvectors from node to all
    r = broadcast(world, rvectors, node)
    # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
    v = np.exp(-1j * np.tensordot(r, gvectors, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    dummy = np.tensordot(v, wfns, ((1),(axis)))
    # reduce across processes
    if node == comm.rank: result = reduce(comm, dummy, lambda x,y: x+y, node)
    else: reduce(comm, dummy, lambda x,y: x+y, node)

  return result

def create_structure():
  from numpy import matrix, array
  from lada.crystal import fill_structure, Structure, Atom

  # creating unrelaxed structure.
  cell = matrix( [ [  0, 0.5, 0.5],
                   [0.5,   0, 0.5],
                   [0.5, 0.5,   0] ] )
  Si = fill_structure(cell)
  for atm in Si.atoms: atm = "Si"

  Ge = fill_structure(cell)
  for atm in Ge.atoms: atm = "Ge"

  return Si, Ge

from sys import exit
from math import ceil, sqrt, pi
from os.path import join, exists
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world, all_reduce, broadcast
from lada.vff import Vff
from lada.escan import Escan, method, nb_valence_states as nbstates, potential, Wavefunctions, to_realspace
from lada.crystal import deform_kpoint
from lada.opt.tempdir import Tempdir
from lada.physics import a0

from numpy import array
import numpy as np
# file with escan and vff parameters.
input = "test_input/input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates unrelaxed structure and  known relaxed structure (for testing).
Si, Ge = create_structure()

# creates vff using parameters in input file. 
# vff will run on all processes (world).
vff = Vff(input, world)

# launch vff and checks output.
Si_relaxed, stress = vff(Si)
Ge_relaxed, stress = vff(Ge)

# creates escan functional using parameters in input file.
# escan will work on all processes.
escan = Escan(input, world)
# sets the FFT mesh
escan.genpot.mesh = (4*20, 20, int(ceil(sqrt(0.5) * 20)))  
escan.genpot.mesh = tuple([int(ceil(sqrt(0.5) * 20)) for i in range(3)])  
escan.genpot.mesh = (16,16,16)
# makes sure we keep output directories.
escan.destroy_directory = False
# calculations are performed using the folded-spectrum method
escan.method = method.folded

# some kpoints
X = array( [0,0,1], dtype="float64" ), "X" # cartesian 2pi/a
G = array( [0,0,0], dtype="float64" ), "Gamma"
L = array( [0.5,0.5,0.5], dtype="float64" ), "L"
W = array( [1, 0.5,0], dtype="float64" ), "W"

# launch pescan for different jobs.
for (kpoint, name), structure, relaxed in [ (G, Si, Si_relaxed) ]:
  with Tempdir(workdir="work", comm=world, keep=True, debug="debug") as tempdir:
    # will save output to directory "name".
    escan.directory = tempdir
    # computes at kpoint of deformed structure.
    escan.kpoint = deform_kpoint(kpoint, structure.cell, relaxed.cell)
    # computing 4 (spin polarized) states.
    escan.nbstates = nbstates(relaxed)
    # divides by two if calculations are not spin polarized.
    if norm(escan.kpoint) < 1e-6 or escan.potential != potential.spinorbit: escan.nbstates /= 2
    # Now just do it.
#   eigenvalues = escan(vff, relaxed)
    eigenvalues = array([-0.91784761, -0.91784761, -0.99130073, 2.42593931])
    # checks expected are as expected. 
    if world.rank == 0: print "Ok - %s: %s -> %s: %s" % (name, kpoint, escan.kpoint, eigenvalues)

    volume = structure.scale / a0("A")
    volume = np.linalg.det(np.matrix(Si_relaxed.cell) * volume)
    with Wavefunctions(escan, [i for i in range(len(eigenvalues))]) as (wfns, gvectors, projs):
      print wfns.shape, gvectors.shape, projs.shape
      rspace, rvectors = to_realspace(wfns, escan.comm)
      print rspace.shape, rvectors.shape
      rvec = broadcast(world, rvectors[1], 0) 
      gvec = broadcast(world, gvectors[1], 0)
      if world.rank == 0: 
        cell = matrix(Si_relaxed.cell* structure.scale /a0("A"))
        kcell = cell.I.T
        print gvec, np.dot(kcell, array([2e0*pi, 0, 0])) 
        print rvec, (cell * matrix([0,0, 1e0/escan.genpot.mesh[0]]).T).T

      for i in range(wfns.shape[1]):
        a = np.matrix(rspace[:,i,0]) # np.multiply(wfns[:,i,j], projs))
        c = all_reduce(world, a*a.T.conjugate(), lambda x,y: x+y)
        a = np.matrix(rspace[:,i,1]) # np.multiply(wfns[:,i,j], projs))
        c += all_reduce(world, a*a.T.conjugate(), lambda x,y: x+y)
        s = all_reduce(world, rspace.shape[0], lambda x,y: x+y)
        if world.rank == 0: print c * volume/s, volume
        u = 0e0 + 0e0*1j
        v = np.exp(-1j * np.dot(gvectors, rvec) )
        u = np.dot(wfns[:,i,0], v)
#       for k in range(wfns.shape[0]):
#         u += np.exp(-1j *np.dot(rvec, gvectors[k, :])) * wfns[k,i,0] # * projs[k]
        u = all_reduce(world, u, lambda x,y:x+y)
        if world.rank == 0: print u, rspace[1,i,0], u.real/rspace[1,i,0].real

      v = np.exp(-1j * np.dot(gvectors, rvec))
      u = np.tensordot(wfns, v, ((0), (0)))
      u = all_reduce(world, u, lambda x,y:x+y)
      if world.rank == 0:
        for test, ref in zip(u.flat, rspace[1,:,:].flat):
          print abs(test - ref) < 1e-12

      r = gtor_fourrier(wfns, rvectors, gvectors, escan.comm)
      assert len(r.flat) == len(rspace.flat)
      for rtest, rtrue in zip(r.flat, rspace.flat):
        assert abs(rtest-rtrue) < 1e-12,\
               RuntimeError("Did not checkout %e != %e" %(rtest, rtrue))
