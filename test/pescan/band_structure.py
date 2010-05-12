import pickle
from sys import exit
from os.path import join
from numpy import matrix, array
from boost.mpi import world
from lada.opt import read_input
from lada.escan import Escan, nb_valence_states, soH
from lada.vff import Vff
from lada.crystal import fill_structure, sort_layers

# reads input file.
global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
input = read_input("input.py", global_dict=global_dict)

# creating unrelaxed structure.
cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
structure = sort_layers(fill_structure(cell), array([0,0,2.5]))
for i, atom in enumerate(structure.atoms):
  atom.type = "Si" if i % 10 < 6 else "Ge"

# some kpoints + associated name
X = array( [0,0,1], dtype="float64" )
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" )
W = array( [1, 0.5,0], dtype="float64" )

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
kpoints = [ (X, G), (G, L) ]
density = 15 / min( norm(X), norm(L), norm(W) )

result = band_structure( input.escan, structure, kpoints, density, 
                         outdir = join("results", "Si"),
                         eref   = None, 
                         nbstates = nb_valence_states(structure) + 4 )
  
if world.rank == 0:
  with open(join("results", "pickle"), "w") as file:
    pickle.dump(result, file) 

