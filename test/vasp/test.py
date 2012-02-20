import lada
lada.lada_with_mpi = False
from lada.crystal import fill_structure
from lada.vasp import read_input
from lada.mpi import NullComm

input = read_input("input.py")

structure = fill_structure(input.lattice.cell, input.lattice)

input.relaxer.vasp.launch_as_library = False
input.relaxer.vasp.program = "vasp-4.6"
print ">>>"
comm = NullComm(4)
print "0. size", comm.size
input.relaxer(structure, outdir="results", comm=comm)
