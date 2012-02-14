import lada
# lada.lada_with_mpi = False
from lada.crystal import fill_structure
from lada.vasp import read_input
from lada.vasp.emass import reciprocal
from lada.mpi import world

input = read_input("input.py")

structure = fill_structure(input.lattice.cell, input.lattice)

input.relaxer.vasp.launch_as_library = True
input.relaxer.vasp.program = "vasp-4.6"
input.relaxer.vasp.vasp_library = "libvasp-4.6.so"
comm = world


result = reciprocal(input.relaxer, structure, outdir="results", comm=comm)
