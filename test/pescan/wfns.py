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
from math import ceil, sqrt
from os.path import join, exists
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world
from lada.vff import Vff
from lada.escan import Escan, method, nb_valence_states as nbstates, potential, wavefunctions
from lada.crystal import deform_kpoint
from lada.opt.tempdir import Tempdir

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
  with Tempdir(workdir="work", keep=True) as tempdir:
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
    # checks expected are as expected. 
    if world.rank == 0: print "Ok - %s: %s -> %s: %s" % (name, kpoint, escan.kpoint, eigenvalues)

    gpoints, wfns = wavefunctions(escan, [i for i in range(len(eigenvalues))])
    print wfns.shape
