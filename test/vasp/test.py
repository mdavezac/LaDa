from shutil import copy as copyfile
from lada.crystal import fill_structure
from lada.opt.changedir import Changedir
from lada.vasp import read_input

input = read_input("input.py")

structure = fill_structure(input.lattice.cell, input.lattice)

input.relaxer(structure, outdir="results")
