# davezac@hopper08:~/SiGe/old/jwluo/SiGe_SL_001/Si6Ge4_on_Ge/vff
from operator import attrgetter
from numpy import array, dot
from numpy.linalg import det
from quantities import angstrom, eV
from boost.mpi import world
from lada.escan import read_input, bandgap
from lada.crystal import fill_structure, sort_layers, nb_valence_states
from lada.crystal.binary import zinc_blende
from lada.physics import a0, reduced_reciprocal_au
from numpy import pi

input = read_input("input.py")
cell = [[2.5,0,0],[0.5,0.5,0.5],[0,-0.5,0.5]]
subs = {'A':'Si', 'B':'Ge'}
structure = zinc_blende().to_structure(cell, subs)
structure.scale = float(5.34594034 / 0.5 * a0.rescale(angstrom))
structure = sort_layers(structure, array([2.5,0,0]))
for i, atom in enumerate(structure.atoms):
  atom.type = 'Si' if i % 10 < 6 else 'Ge'

input.escan.fft_mesh = 50, 14, 45
input.escan.vff.direction = None

result = bandgap( input.escan, structure, comm=world, eref=None, 
                  outdir="result/dipoles", references = (-0.2, 0.2),
                  nbstates=6, direction = array([1., 0,0]) )
volume = det(result.structure.cell*result.structure.scale/a0.rescale(angstrom))
cbms = result.extract_cbm.gwfns 
# cbms = sorted(cbms, key=attrgetter("eigenvalue"))
vbms = result.extract_vbm.gwfns 
# vbms = sorted(vbms, key=attrgetter("eigenvalue"))
units = pi ** 2 # / 4. #float(a0.rescale(angstrom) * a0.rescale(angstrom))
dips = []
for e0, e1, u in result.dipole():
  dips.append((e0, e1, (u * u.conjugate()).real))
check = [(array(0.168325) * eV, array(-0.49477701000000002) * eV, array([  1.43562357e-04,   1.04201091e-05,   1.95678721e-05]) * 1/a0**2), (array(0.168325) * eV, array(-0.49477701000000002) * eV, array([  4.20854765e-04,   1.67207771e-05,   7.57293883e-06]) * 1/a0**2), (array(0.168325) * eV, array(-0.49477701000000002) * eV, array([  4.20854765e-04,   1.67207771e-05,   7.57293883e-06]) * 1/a0**2), (array(0.168325) * eV, array(-0.49477701000000002) * eV, array([  1.43562357e-04,   1.04201091e-05,   1.95678721e-05]) * 1/a0**2)]
for (e0, e1, a), (E0, E1, A) in zip(dips, check):
  assert abs(e0-E0) < 1e-6
  assert abs(e1-E1) < 1e-6
  assert all(abs(a-A) < 1e-6)
