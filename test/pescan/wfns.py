
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
  """ g-space to r-space fourrier transform of wavefunctions.
  
      @param wavefunctions: an numpy array of wavefunctions.
      @param rvectors: a two-dimensional array of r-space vectors, with each
        row a position. The return r-space wavefunctions will be given with
        respect to these points. Each process should have different r-space
        vectors. Otherwise use another implementation.
      @param gvectors: a two-dimensional array of g-space vectors, with each
        row a (g-)position. The input wavefunctions should be given with
        respect to these points, in the same order, etc.
      @param comm: communicator over which the wavefunctions are distributed.
        The return wavefunctions will also be dirstributed over these
        processes.
      @params axis: axis over which the wavefunctions are deployed, eg axis
        of L{wavefunctions} which corresponds to L{gvectors}. Independent
        (though simultaneous, implementation wise) fourrier transform will be
        performed over this axis for all other axis. 0 by default.
      @type axis: integer
  """
  import numpy as np
  from boost.mpi import broadcast, reduce

  result = None
  for node in range(comm.size):
    # sends rvectors from node to all
    r = broadcast(world, rvectors, node)
    # computes all exponentials exp(-i r.g), with r in first dim, and g in second.
    v = np.exp(-1j * np.tensordot(r, gvectors, ((1),(1))))
    # computes fourrier transform for all wavefunctions simultaneously.
    dummy = np.tensordot(v, wavefunctions, ((1),(axis)))
    # reduce across processes
    if node == comm.rank: result = reduce(comm, dummy, lambda x,y: x+y, node).copy()
    else: reduce(comm, dummy, lambda x,y: x+y, node)

  return result
def rtog_fourrier(wavefunctions, rvectors, gvectors, comm, axis=0):
  """ r-space to g-space fourrier transform of wavefunctions.
  
      @param wavefunctions: an numpy array of wavefunctions.
      @param rvectors: a two-dimensional array of r-space vectors, with each
        row a position. The return r-space wavefunctions will be given with
        respect to these points. Each process should have different r-space
        vectors. Otherwise use another implementation.
      @param gvectors: a two-dimensional array of g-space vectors, with each
        row a (g-)position. The input wavefunctions should be given with
        respect to these points, in the same order, etc.
      @param comm: communicator over which the wavefunctions are distributed.
        The return wavefunctions will also be dirstributed over these
        processes.
      @params axis: axis over which the wavefunctions are deployed, eg axis
        of L{wavefunctions} which corresponds to L{gvectors}. Independent
        (though simultaneous, implementation wise) fourrier transform will be
        performed over this axis for all other axis. 0 by default.
      @type axis: integer
  """
  import numpy as np
  from boost.mpi import broadcast, reduce, all_reduce

  result = None
  for node in range(comm.size):
    # sends rvectors from node to all
    g = broadcast(world, gvectors, node)
    # computes all exponentials exp(-i r.g), with g in first dim, and r in second.
    v = np.exp(1j * np.tensordot(g, rvectors, ((1),(1))))
    # gets normalization factor.
    norm = all_reduce(comm, rvectors.shape[0], lambda x,y:x+y)
    # computes fourrier transform for all wavefunctions simultaneously.
    dummy = np.tensordot(v, wavefunctions, ((1),(axis))) / float(norm)
    # reduce across processes
    if node == comm.rank: result = reduce(comm, dummy, lambda x,y: x+y, node).copy()
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
for (kpoint, name), structure, relaxed in [(G, Si, Si_relaxed), (G, Ge, Ge_relaxed) ]:
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
    eigenvalues = escan(vff, relaxed)
#   eigenvalues = array([-0.94213329, -0.94213329, -0.94213329, 2.46059562])
    # checks expected are as expected. 
    if world.rank == 0: print "Ok - %s: %s -> %s: %s" % (name, kpoint, escan.kpoint, eigenvalues)

    volume = structure.scale / a0("A")
    volume = np.linalg.det(np.matrix(Si_relaxed.cell) * volume)
    with Wavefunctions(escan, [i for i in range(len(eigenvalues))]) as (wfns, gvectors, projs):
      rspace, rvectors = to_realspace(wfns, escan.comm)

      other = rtog_fourrier(rspace, rvectors, gvectors, escan.comm)
      assert len(other.flat) == len(wfns.flat)
      for gtest, gtrue in zip(other.flat, wfns.flat):
        assert abs(gtest-gtrue) < 1e-12,\
               RuntimeError("Did not checkout %e != %e" %(gtest, gtrue))
      print "rtog passed", escan.comm.rank

      other = gtor_fourrier(wfns, rvectors, gvectors, escan.comm)
      assert len(other.flat) == len(rspace.flat)
      for rtest, rtrue in zip(other.flat, rspace.flat):
        assert abs(rtest-rtrue) < 1e-12,\
               RuntimeError("Did not checkout %e != %e" %(rtest, rtrue))
      print "gtor passed", escan.comm.rank

