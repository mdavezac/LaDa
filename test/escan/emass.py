from math import fabs as abs
from numpy import array
from lada.crystal import fill_structure
from lada.escan import read_input, extract_bg
from lada.escan.emass import Functional

input = read_input("input.py", namespace = {"Escan": Functional})

G = array([0,0,0], dtype="float64")
X = array([0,0,1], dtype="float64")
structure = fill_structure(input.vff.lattice.cell, input.vff.lattice)
input.escan.nbstates = len(structure.atoms) * 4 + 4
input.escan.tolerance = 1e-12

# result = input.escan( structure, direction=(0,0,1), outdir="results/emass/100", \
#                       do_relax_kpoint=False, type="e" )
# print result.mass
# assert abs(result.mass[0] - 0.4381) < 0.01 # Gamma conduction emass 

result = input.escan( structure, direction=(0,0,1), outdir="results/hmass/100", \
                      do_relax_kpoint=False, type="h", bandgap=extract_bg("results/emass/100") )
print result.mass
assert abs(result.mass[0] - 0.2769) < 0.01 # Gamma conduction heavy hmass (100 direction)
assert abs(result.mass[2] - 0.2059) < 0.01 # Gamma conduction light hmass (100 direction)

result = input.escan( structure, direction=(1,1,1), outdir="results/hmass/111", \
                      do_relax_kpoint=False, type = "h", bandgap=result.extract_bg )
print result.mass
assert abs(result.mass[0] - 0.6885) < 0.01 # Gamma conduction heavy hmass (111 direction)
assert abs(result.mass[2] - 0.1460) < 0.01 # Gamma conduction light hmass (111 direction)

# result = reciprocal( input.escan, structure, X, outdir="results/emass/100", \
#                      order = 2, do_relax_kpoint=False )
# assert abs(1e0 / result[0][2,8] - 0.4381) < 0.01 # Gamma conduction emass 
# assert abs(-1e0 / result[0][2,7] - 0.2769) < 0.01 # Gamma conduction heavy hmass (100 direction)
# assert abs(-1e0 / result[0][2,5] - 0.2059) < 0.01 # Gamma conduction light hmass (100 direction)
# result = reciprocal( input.escan, structure, array([1e0, 1, 1]), outdir="results/emass/111", \
#                      order = 2, do_relax_kpoint=False )
# assert abs(-1e0 / result[0][2,7] - 0.6885) < 0.01 # Gamma conduction heavy hmass (111 direction)
# assert abs(-1e0 / result[0][2,5] - 0.1460) < 0.01 # Gamma conduction light hmass (111 direction)
