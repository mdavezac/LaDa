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
  Si = fill_structure(lattice.cell)
  Ge = fill_structure(lattice.cell)
  for atom in Si.atoms: atom.type = "Si"
  for atom in Ge.atoms: atom.type = "Ge"
  return Si, Ge

import cPickle
from sys import exit
from math import ceil, sqrt
from os.path import join, exists
from os import getcwd
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world
from lada.vff import Vff
from lada.escan import Escan, method, nb_valence_states as nbstates, potential, derivatives
from lada.crystal import deform_kpoint
from lada.opt.tempdir import Tempdir
from lada.crystal import Atoms, Atom
from sys import exit

# file with escan and vff parameters.
input = "test_input/input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates unrelaxed structure and  known relaxed structure (for testing).
Si, Ge = create_structure()

# creates vff using parameters in input file. 
# vff will run on all processes (world).
vff = Vff(input, world)

# launch vff 
Si_relaxed, stress = vff(Si)
Ge_relaxed, stress = vff(Ge)

# creates escan functional using parameters in input file.
# escan will work on all processes.
escan = Escan(input, world)
# sets the FFT mesh
escan.genpot.mesh = tuple(int(ceil(sqrt(0.5) * 20)) for u in range(3))  

# some kpoints + associated emass direction.
X = array( [0,0,1], dtype="float64" ), "X"
G = array( [0,0,0], dtype="float64" ), "Gamma"
L = array( [0.5,0.5,0.5], dtype="float64" ) "L"
W = array( [1, 0.5,0], dtype="float64" ), "W"

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
jobs = [ (G, X), (X, L), (L, G), (G, W) ]
density = 10 / norm(X[0]) 

Si_bandstructure = BandStructure( Si_relaxed, escan, vff, jobs, density, 
                                  directory = join("work", "Si"),
                                  method = method.full_diagonalization,
                                  destroy_directory = False )
Ge_bandstructure = BandStructure( Ge_relaxed, escan, vff, jobs, density, 
                                  directory = join("work", "Si"),
                                  method = method.full_diagonalization,
                                  destroy_directory = False )


with open(join("work", "pickled_sige"), "r") as file:
  cPickle.dump((Si_bandstructure, Ge_bandstructure), file) 
