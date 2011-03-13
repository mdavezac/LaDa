# davezac@hopper08:~/SiGe/old/jwluo/SiGe_SL_001/Si6Ge4_on_Ge/vff
from operator import attrgetter
from numpy import array, dot
from numpy.linalg import det
from quantities import angstrom, eV
from lada.escan import read_input, bandgap
from lada.crystal import fill_structure, sort_layers, nb_valence_states, Structure
from lada.crystal.binary import zinc_blende
from lada.physics import a0, reduced_reciprocal_au
from lada.mpi import world
from numpy import pi, multiply

input = read_input("input.py")
cell = [[2.5,0,0],[0.5,0.5,0.5],[0,-0.5,0.5]]
subs = {'A':'Si', 'B':'Ge'}

structure = zinc_blende().to_structure(cell, subs)
structure.scale = float(5.33847592 / 0.5 * a0.rescale(angstrom))
structure = sort_layers(structure, array([2.5,0,0]))
for i, atom in enumerate(structure.atoms):
  atom.type = 'Si' if i % 10 < 6 else 'Ge'

result = bandgap( input.escan, structure, eref=None, 
                  outdir="results/dipoles", references = (-0.2, 0.2),
                  nbstates=4, direction = array([1., 0,0]),
                  fft_mesh = (50,14,14) )

tot = array([0e0] * 3)  / a0 / a0
for e0, e1, u in result.dipole(degeneracy = 5e1 * input.escan.tolerance, attenuate=False):
  tot += (u * u.conjugate()).real
if result.extract_vbm.nnodes == 1: # different from MPI! escan bug...
  check = array([1.13781377e-03,   5.90701883e-05,   5.90701874e-05]) * 1/a0**2
else: 
  check = array([ 0.00033923,  0.00153695,  0.00153694]) * 1/a0**2
assert all( abs(check - tot) < min(abs(check))*1e-3 ), abs(check-tot)
