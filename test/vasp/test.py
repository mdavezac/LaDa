from shutil import copy as copyfile
from boost.mpi import world
from lada.opt import read_input
from lada.crystal import fill_structure
from lada.opt.changedir import Changedir
from lada.vasp import Vasp, specie
from lada.vasp.methods import RelaxCellShape

input_dict = { "Vasp": Vasp, "U": specie.U, "nlep": specie.nlep }
input = read_input("input.py", input_dict)

structure = fill_structure(input.lattice.cell, input.lattice)

input.vasp(structure, comm=world, outdir="results/0")
with Changedir("results/1") as pwd: pass
copyfile("results/0/CONTCAR", "results/1/CONTCAR")
input.vasp(structure, comm=world, outdir="results/1")

