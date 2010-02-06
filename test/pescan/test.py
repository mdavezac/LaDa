def create_zb_lattice():
  from numpy import array as np_array
  from lada.crystal import Lattice, Site

  lattice = Lattice()
  lattice.cell = np_array( [ [  0, 0.5, 0.5],
                             [0.5,   0, 0.5],
                             [0.5, 0.5,   0] ] )
  lattice.sites.append( Site( np_array([0.0,0,0]), ["Si", "Ge"] ) )
  lattice.sites.append( Site( np_array([0.25,0.25,0.25]), ["Si", "Ge"] ) )
  lattice.scale = 5.45
  lattice.find_space_group()
  lattice.set_as_crystal_lattice() 
  return lattice

def create_structure():
  from numpy import matrix, array
  from lada.crystal import fill_structure, Structure, Atom

  # creating unrelaxed structure.
  cell = matrix( [ [4,  0, 0.5],
                   [0,  1,   0],
                   [0,  0, 0.5] ] )
  structure = fill_structure(cell)
  N = len(structure.atoms)
  for i in range(N/2):
    site = structure.atoms[i].site
    structure.atoms[i].type = structure.lattice.sites[site].type[0]
  for i in range(N/2, N):
    site = structure.atoms[i].site
    structure.atoms[i].type = structure.lattice.sites[site].type[1]

  # The relaxed structure should look like this.
  result = Structure()
  result.cell = matrix( [ [ 4.07074, -0.00000,  0.50814],  
                          [-0.00000,  1.01448,  0.00000],  
                          [-0.00558,  0.00000,  0.50814] ] ) 
  result.atoms.append( Atom(array([0.00000,  0.00000,  0.00000]), "Si", 0) )
  result.atoms.append( Atom(array([0.25407,  0.24762,  0.25407]), "Si", 1) )
  result.atoms.append( Atom(array([3.56953,  0.50575, -0.01252]), "Si", 0) )
  result.atoms.append( Atom(array([3.82360,  0.75634,  0.24155]), "Si", 1) )
  result.atoms.append( Atom(array([3.06705, -0.00149, -0.01818]), "Si", 0) )
  result.atoms.append( Atom(array([3.32112,  0.24911,  0.23589]), "Si", 1) )
  result.atoms.append( Atom(array([2.56584,  0.50724, -0.02511]), "Si", 0) )
  result.atoms.append( Atom(array([2.81991,  0.75486,  0.22896]), "Si", 1) )
  result.atoms.append( Atom(array([2.05764, -0.00658, -0.02507]), "Ge", 0) )
  result.atoms.append( Atom(array([2.31172,  0.25420,  0.22901]), "Ge", 1) )
  result.atoms.append( Atom(array([1.54094,  0.50155, -0.01651]), "Ge", 0) )
  result.atoms.append( Atom(array([1.79501,  0.76054,  0.23756]), "Ge", 1) )
  result.atoms.append( Atom(array([1.02490, -0.00569, -0.00861]), "Ge", 0) )
  result.atoms.append( Atom(array([1.27897,  0.25330,  0.24547]), "Ge", 1) )
  result.atoms.append( Atom(array([0.50819,  0.50066, -0.00005]), "Ge", 0) )
  result.atoms.append( Atom(array([0.76227,  0.76143,  0.25402]), "Ge", 1) )
  result.energy = 0.0927878437055998

  return structure, result

def deformed_kpoint(ocell, dcell, kpoint):
  from numpy import dot, matrix
  k = dot(ocell.T, kpoint)
  for i in range(3):
    k[i] = float(k[i]) - float( int(k[i]) )
    if k[i] < 0e0: k[i] += 1e0
    if k[i] > 1e0-1e-6: k[i] = 0e0
  return matrix(dcell.T).I * matrix(k).T
  

from sys import exit
from math import ceil, sqrt
from os.path import join, exists
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world
from lada.vff import Vff
from lada.escan import Escan, method, nb_valence_states as nbstates, potential

# file with escan and vff parameters.
input = "input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates unrelaxed structure and  known relaxed structure (for testing).
structure, result = create_structure()

# creates vff using parameters in input file. 
# vff will run on all processes (world).
vff = Vff(input, world)

# launch vff and checks output.
relaxed, stress = vff(structure)
if world.rank == 0:
  diff = 0e0
  for a, b in zip(relaxed.atoms, result.atoms):
    assert a.type == b.type
    assert a.site == b.site
    diff += dot(a.pos-b.pos, a.pos-b.pos)

  assert diff / float(len(structure.atoms)) < 1e-8, diff 
  assert abs(relaxed.energy - result.energy) < 1e-8, abs(relaxed.energy - result.energy) 

# creates escan functional using parameters in input file.
# escan will work on all processes.
escan = Escan(input, world)
# sets the FFT mesh
escan.genpot.mesh = (4*20, 20, int(ceil(sqrt(0.5) * 20)))  
# makes sure we keep output directories.
escan.destroy_directory = False
# calculations are performed using the folded-spectrum method
escan.method = method.folded

# some kpoints
X = array( [0,0,1], dtype="float64" ) # cartesian 2pi/a
G = array( [0,0,0], dtype="float64" )
L = array( [0.5,0.5,0.5], dtype="float64" ) 
W0 = array( [1, 0.5,0], dtype="float64" ) 
W1 = array( [1, 0,0.5], dtype="float64" ) 
W2 = array( [0, 1,0.5], dtype="float64" ) 

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
jobs = [\
         (G,   "VBM", -0.4, array([-0.47992312, -0.67148097])), # at gamma, code uses Krammer degeneracy
         (G, "Gamma",  0.4, array([ 0.47368306,  0.49199994])), # at gamma, code uses Krammer degeneracy
         (X,     "X",  0.4, array([ 0.51468608,  0.51479076, 0.5148467 , 0.5149207 ])),
         (L,     "L",  0.4, array([ 0.72789198,  0.72789198, 0.73165765, 0.73165765])),
         (W1,   "W1",  0.4, array([ 0.89170814,  0.89170822, 0.96097565, 0.96097601])),
         (W2,   "W2",  0.4, array([ 0.89174454,  0.89174462, 0.9608853 , 0.96088566]))
       ]
# launch pescan for different jobs.
for kpoint, name, ref, expected_eigs in jobs:
  # will save output to directory "name".
  escan.directory = name
  # computes at kpoint of deformed structure.
  escan.kpoint = deformed_kpoint(structure.cell, relaxed.cell, kpoint)
  # computing 4 (spin polarized) states.
  escan.nbstates = 4 # + nbtates(relaxed)  would be all valence states + 4
  # divides by two if calculations are not spin polarized.
  if norm(escan.kpoint) < 1e-6 or escan.potential != potential.spinorbit: escan.nbstates /= 2
  # sets folded method's energy reference
  escan.reference = ref
  # Now just do it.
  eigenvalues = escan(vff, relaxed)
  # checks expected are as expected. 
  assert norm( eigenvalues - expected_eigs ) < 1e-6, "%s\n%s" % (eigenvalues, expected_eigs)
  # And print.
  if world.rank == 0: print "Ok - %s: %s -> %s: %s" % (name, kpoint, escan.kpoint, eigenvalues)
  
