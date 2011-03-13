from shutil import copy as copyfile
from lada.crystal import fill_structure
from lada.opt.changedir import Changedir
from lada.vasp import read_input

input = read_input("input.py")

structure = fill_structure(input.lattice.cell, input.lattice)

input.vasp(structure, outdir="results/0")
with Changedir("results/1") as pwd: pass
copyfile("results/0/CONTCAR", "results/1/CONTCAR")
input.vasp(structure, outdir="results/1")

