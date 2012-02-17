from numpy import dot
from lada.crystal import fill_structure, FreezeAtom
from lada.vasp import read_input
from lada.vasp.emass import reciprocal
from lada.mpi import world

input = read_input("input.py")

structure = input.lattice.to_structure() 
for atom in structure.atoms: atom.type = "Ge"
structure.scale = 5.65
# structure = fill_structure(dot(input.lattice.cell, [[-1,0,0],[1,0,1],[1,1,0]]), input.lattice)
# for atom in structure.atoms: atom.freeze = FreezeAtom.x

input.relaxer.vasp.ispin = 1
input.relaxer.vasp.kpoints    = "Automatic generation\n0\nGamma\n4 4 4\n0 0 0"
input.relaxer.vasp.launch_as_library = True
input.relaxer.vasp.program = "vasp-4.6"
input.relaxer.vasp.vasp_library = "libvasp-4.6.so"
comm = world


result = reciprocal(input.relaxer, structure, outdir="results", stepsize=0.1, order=4, comm=comm)
