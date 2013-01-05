# davezac@hopper08:~/SiGe/old/jwluo/SiGe_SL_001/Si6Ge4_on_Ge/vff
from numpy import array
from quantities import angstrom
from pylada.escan import read_input, bandgap
from pylada.crystal import sort_layers
from pylada.crystal.binary import zinc_blende
from pylada.physics import a0

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
check = array([1.13753549e-03,   5.87932303e-05,   5.87932303e-05]) * 1/a0**2
for a, b in zip(abs(check - tot), abs(check)*1e-2): assert a < max(b, 1e-6), (a, b)
