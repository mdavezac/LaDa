# davezac@hopper08:~/SiGe/old/jwluo/SiGe_SL_001/Si6Ge4_on_Ge/vff
from operator import attrgetter
from numpy import array, dot
from numpy.linalg import det
from quantities import angstrom
from boost.mpi import world
from lada.escan import read_input, bandgap
from lada.crystal import fill_structure, sort_layers, nb_valence_states
from lada.crystal.binary import zinc_blende
from lada.physics import a0

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
print result.extract_vbm.eigenvalues
print result.extract_cbm.eigenvalues
volume = det(result.structure.cell*result.structure.scale/a0.rescale(angstrom))
print "volume real space:", volume
cbms = result.extract_cbm.gwfns 
cbms = sorted(cbms, key=attrgetter("eigenvalue"))
vbms = result.extract_vbm.gwfns 
vbms = sorted(vbms, key=attrgetter("eigenvalue"))
for vbm in vbms[::2]:
  for cbm in cbms:
    dip = vbm.braket(result.extract_vbm.gvectors[:,0], cbm, attenuate=False)
    print float((dip * dip.conjugate()).real) * volume, 
  print
# if world.rank == 0:
#   dips = result.dipole()
#   for dip in dips:
#     print dip[2]
