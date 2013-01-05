from numpy import dot
from pylada.crystal import fill_structure, FreezeAtom
from pylada.vasp import read_input
from pylada.vasp.emass import reciprocal
from pylada.mpi import world

input = read_input("input.py")

structure = fill_structure(dot(input.lattice.cell, [[-1,0,0],[1,0,1],[1,1,0]]), input.lattice)
for atom in structure.atoms: atom.freeze = FreezeAtom.x

input.relaxer.vasp.launch_as_library = True
input.relaxer.vasp.program = "vasp-4.6"
input.relaxer.vasp.vasp_library = "libvasp-4.6.so"
comm = world


result = reciprocal(input.relaxer, structure, outdir="results", comm=comm)
