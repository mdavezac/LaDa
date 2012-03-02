import lada
lada.lada_with_mpi = False
from lada.vasp import read_input

input = read_input("input.py")
input.structure.name = "has a name"

input.vasp(input.structure, outdir="results", comm={'n': 2, 'ppn': 1})
