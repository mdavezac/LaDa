""" Computes dipole matrix elements and oscillator strengths. """
import pickle
from sys import exit
from os import getcwd
from os.path import join
from numpy import matrix, array, conjugate, transpose, multiply, dot
from numpy.linalg import norm, det
from boost.mpi import world
from lada.opt import read_input, redirect
from lada.opt.changedir import Changedir
from lada.escan import Escan, soH, bandgap, Extract
from lada.escan._bandgap import ExtractRefs
from lada.escan._wfns import rtog_fourrier
from lada.vff import Vff
from lada.crystal import fill_structure, sort_layers, FreezeCell, nb_valence_states 
from lada.physics import a0

# reads input file.
global_dict={"Vff": Vff, "Escan": Escan, "nb_valence_states": nb_valence_states, "soH": soH}
input = read_input("input.py", global_dict=global_dict)

# creating unrelaxed structure.
input.vff.lattice.set_as_crystal_lattice()
cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
structure = sort_layers(fill_structure(cell), array([0,0,2.5]))
for i, atom in enumerate(structure.atoms):
  atom.type = "Si" if i % 10 < 6 else "Ge"
# structure = fill_structure(input.escan.vff.lattice.cell)
# structure.atoms[0].type = "Si"
structure.scale = 5.65
structure.freeze = FreezeCell.a0 | FreezeCell.a1
input.escan.fft_mesh  = 14, 14, 50

out = bandgap( input.escan, structure,\
               outdir=join("results", "osc"),\
               references=(0.1, -0.4),\
               nbstates = 4,\
               comm=world )

# gvectors = out.extract_cbm.gvectors
# vol = det(out.extract_vbm.structure.cell * out.extract_vbm.structure.scale / a0("A"))
# for i, awfn in enumerate(out.extract_vbm.gwfns):
#   for j, bwfn in enumerate(out.extract_cbm.gwfns):
#     a = awfn.braket(transpose(gvectors), bwfn, attenuate=True) 
#     if world.rank == 0:
#       print multiply(a, a.conjugate()) * vol * vol, dot(a, a.conjugate()).real * (vol * vol)
#   break
#   print 
osc = out.oscillator_strength() 
if world.rank == 0: print osc

# results =   dipole_matrix_elements(out.extract_vbm, out.extract_cbm) \
#           * det(out.extract_vbm.structure.cell * out.extract_vbm.structure.scale / a0("A"))
# print  results
