import pickle
from sys import exit
from os.path import join
from numpy import matrix, array
from numpy.linalg import norm
from boost.mpi import world
from lada.opt import read_input
from lada.escan import Escan, soH, bandgap
from lada.vff import Vff
from lada.crystal import fill_structure, sort_layers, FreezeCell, nb_valence_states

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

out = bandgap( input.escan, structure,\
               outdir=join("results", "BG"),\
               comm=world, references = (0.2, -0.5 ) )
print "%f - %f = %f " % (out.cbm, out.vbm, out.bandgap)
