import pickle
from sys import exit
from os.path import join
from numpy import matrix, array
from numpy.linalg import norm
from boost.mpi import world
from lada.opt import read_input
from lada.escan import Escan, nb_valence_states, soH, band_structure
from lada.vff import Vff
from lada.crystal import fill_structure, sort_layers, FreezeCell

# reads input file.
global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
input = read_input("input.py", global_dict=global_dict)

# creating unrelaxed structure.
input.vff.lattice.set_as_crystal_lattice()
cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
structure = sort_layers(fill_structure(cell), array([0,0,2.5]))
for i, atom in enumerate(structure.atoms):
  atom.type = "Si" if i % 10 < 6 else "Ge"
structure.scale = 5.65
structure.freeze = FreezeCell.a0 | FreezeCell.a1
input.escan.fft_mesh  = 14, 14, 50

structure = fill_structure(input.vff.lattice.cell)
structure.atoms[0].type = "Si"
input.escan.fft_mesh  = 14, 14, 14
structure.scale = 5.45

# some kpoints + associated name
X = array( [1,0,0], dtype="float64" )
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" )
W = array( [0, 0.5,1], dtype="float64" )

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
kpoints = [ (X, G), (G, L) ]
density = 20 / min( norm(X), norm(L), norm(W) )

result = band_structure( input.escan, structure, kpoints, density, 
                         outdir = "results",
                         eref   = None, 
                         nbstates = nb_valence_states(structure) + 4,
                         pools = 4)
  
if world.rank == 0:
  with open(join("results", "pickle"), "w") as file:
    pickle.dump(result, file) 

