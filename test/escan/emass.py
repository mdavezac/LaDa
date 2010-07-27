from numpy import array
from boost.mpi import world
from lada.crystal import fill_structure
from lada.vff import Vff
from lada.escan import Escan, soH, nb_valence_states
from lada.escan.derivatives import reciprocal
from lada.opt import read_input

input = read_input("input.py", local_dict = {"Vff": Vff, "Escan": Escan, "soH": soH})

G = array([0,0,0], dtype="float64")
X = array([0,0,1], dtype="float64")
structure = fill_structure(input.vff.lattice.cell, input.vff.lattice)
input.escan.nbstates = nb_valence_states(structure) + 4
input.escan.tolerance = 1e-12

print nb_valence_states(structure) 

result = reciprocal( input.escan, structure, X, outdir="work/emass/100", \
                     comm = world, order = 2, _dont_deform_kpoint=True )
print result
result = reciprocal( input.escan, structure, array([1e0, 1, 1]), outdir="work/emass/111", \
                     comm = world, order = 2, _dont_deform_kpoint=True )
print result
