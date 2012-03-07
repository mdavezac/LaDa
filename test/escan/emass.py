import sys
import DLFCN
sys.setdlopenflags(DLFCN.RTLD_NOW)
from math import fabs as abs
from numpy import array
from lada.crystal import fill_structure
from lada.escan import read_input
from lada.escan.emass import Functional

input = read_input("input.py", namespace = {"Escan": Functional})

G = array([0,0,0], dtype="float64")
X = array([0,0,1], dtype="float64")
structure = fill_structure(input.vff.lattice.cell, input.vff.lattice)
input.escan.nbstates = len(structure.atoms) * 4 + 4
input.escan.tolerance = 1e-12

orig = input.escan( structure, direction=(0,0,1), outdir="results/emass/100", \
                      do_relax_kpoint=False, type="e" )
assert abs(orig.mass[0] - 0.4381) < 0.01 # Gamma conduction emass 

result = input.escan( structure, direction=((0,0,1), (1,1,1)), outdir="results/hmass", \
                      do_relax_kpoint=False, type="h", bandgap=orig.extract_bg )
assert abs(result.mass[0][0] - 0.2769) < 0.01 # Gamma conduction heavy hmass (100 direction)
assert abs(result.mass[0][2] - 0.2059) < 0.01 # Gamma conduction light hmass (100 direction)
assert abs(result.mass[1][0] - 0.6885) < 0.01 # Gamma conduction heavy hmass (111 direction)
assert abs(result.mass[1][2] - 0.1460) < 0.01 # Gamma conduction light hmass (111 direction)
