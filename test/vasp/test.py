import lada
lada.lada_with_mpi = False
from shutil import copy as copyfile
from lada.crystal import fill_structure
from lada.opt.changedir import Changedir
from lada.vasp import read_input

input = read_input("input.py")

structure = fill_structure(input.lattice.cell, input.lattice)

input.relaxer.vasp.launch_as_library = False
input.relaxer.vasp.vasp_program = "vasp"
input.relaxer(structure, outdir="results")
