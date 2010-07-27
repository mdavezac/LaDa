from numpy import array
from boost.mpi import world
from lada.crystal import fill_structure
from lada.vff import Vff
from lada.escan import Escan, soH
from lada.escan.derivatives import reciprocal
from lada.opt import read_input

input = read_input("input.py", local_dict = {"Vff": Vff, "Escan": Escan})

G = array([0,0,0], dtype="float64")
X = array([0,0,1], dtype="float64")
structure = fill_structure(input.vff.lattice.cell, input.vff.lattice)
input.escan.nbstates = nb_valence_states(structure) + 8

result = reciprocal( input.escan, structure, X, outdir="work/emass", \
                     comm = world, order = 2 )
print resutl
