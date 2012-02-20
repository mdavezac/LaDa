""" Computes dipole matrix elements and oscillator strengths. """
from os.path import join
from numpy import matrix, array
from boost.mpi import world
from lada.crystal import sort_layers, FreezeCell
from lada.escan import bandgap, read_input

# reads input file.
input = read_input("input.py")

# creating unrelaxed structure.
input.vff.lattice.set_as_crystal_lattice()
cell = matrix([[0.5,-0.5,0],[0.5,0.5,0.5],[0,0,2.5]])
structure = sort_layers(input.lattice.to_structure(cell), array([0,0,2.5]))
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
#     if world.rank == 0:, FreezeCell
#       print multiply(a, a.conjugate()) * vol * vol, dot(a, a.conjugate()).real * (vol * vol)
#   break
#   print 
osc = out.oscillator_strength() 
if world.rank == 0: print osc

# results =   dipole_matrix_elements(out.extract_vbm, out.extract_cbm) \
#           * det(out.extract_vbm.structure.cell * out.extract_vbm.structure.scale / a0("A"))
# print  results
