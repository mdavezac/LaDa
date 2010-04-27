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

def check_emass( structure, escan, vff, direction, order = 1, \
                 nbpoints = None, stepsize = 1e-2, **kwargs ):
  from os import makedirs
  from copy import deepcopy
  from numpy import array
  from lada.escan import eMass
  from lada.opt.changedir import Changedir

  comm = escan.comm
  escan = deepcopy(escan)
  escan.comm = comm
  popthese = []
  for key in kwargs:
    if not hasattr(escan, key): continue
    setattr(escan, key, kwargs[key])
    popthese.append(key)
  for key in popthese: del kwargs[key]

  emass  = eMass()
  emass.order = order
  if nbpoints == None: nbpoints = order+1
  emass.nbpoints = nbpoints + 1
  emass.stepsize = stepsize
  emass.direction = array(direction, dtype="float64")
  emass.nbstates = escan.nbstates
  print emass.nbstates 
  print escan.nbstates 
  emass.kpoint = escan.kpoint
  cell = structure.cell.copy()

  escan.scale = structure
  escan.comm.barrier()
  original = deepcopy(escan.vff_inputfile)
  escan.vff_inputfile = "%s.%i" % (original, escan.comm.rank)
  vff.print_escan_input( escan.vff_inputfile, structure )
  result = emass(escan, cell, structure, escan.reference)
  escan.vff_inputfile = original
  print "s: ", [r[0] for r in result]
  return [ r[1] for r in result ], [r[0] for r in result]
  


from sys import exit
from math import ceil, sqrt
from os.path import join, exists
from os import getcwd
from numpy import dot, array, matrix
from numpy.linalg import norm
from boost.mpi import world
from lada.vff import Vff
from lada.escan import Escan, method, nb_valence_states as nbstates, potential, derivatives
from lada.crystal import deform_kpoint
from lada.opt.tempdir import Tempdir
from lada.crystal import Atoms, Atom
from sys import exit

# file with escan and vff parameters.
input = "test_input/input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates unrelaxed structure and  known relaxed structure (for testing).
structure, result = create_structure()
structure.cell = lattice.cell
structure.atoms = Atoms( [Atom(lattice.sites[0].pos, lattice.sites[0].type[0]),\
                          Atom(lattice.sites[1].pos, lattice.sites[1].type[0])] )

# creates vff using parameters in input file. 
# vff will run on all processes (world).
vff = Vff(input, world)

# launch vff and checks output.
relaxed, stress = vff(structure)
if world.rank == -1:
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

# some kpoints + associated emass direction.
X = array( [0,0,1], dtype="float64" ),       (1,0,0)   
G = array( [0,0,0], dtype="float64" ),       (0,0,1)
L = array( [0.5,0.5,0.5], dtype="float64" ), (1,0,1)
W0 = array( [1, 0.5,0], dtype="float64" ),   (2,0,1) 
W1 = array( [1, 0,0.5], dtype="float64" ),   (2,0,1) 
W2 = array( [0, 1,0.5], dtype="float64" ),   (2,0,1) 

# Each job is performed for a given kpoint (first argument), at a given
# reference energy (third argument). Results are stored in a specific directory
# (second arguement). The expected eigenvalues are given in the fourth argument.
jobs = [\
         # at gamma, code uses Krammer degeneracy
         (X,     "X",   0.4, 8, array([ 0.51468608,  0.51479076, 0.5148467 , 0.5149207 ])),
         (G,   "VBM",  -0.4, 2, array([-0.47992312, -0.67148097])), 
         (G, "Gamma",   0.4, 4, array([ 0.47368306,  0.49199994])), 
         (L,     "L",   0.4, 4, array([ 0.72789198,  0.72789198, 0.73165765, 0.73165765])),
         (W1,   "W1",   0.4, 4, array([ 0.89170814,  0.89170822, 0.96097565, 0.96097601])),
         (W2,   "W2",   0.4, 4, array([ 0.89174454,  0.89174462, 0.9608853 , 0.96088566]))
       ]
# launch pescan for different jobs.
for (kpoint, direction), name, ref, nbstates, expected in jobs:
  # Just do it!
  for i, callme in enumerate([derivatives.reciprocal]): #, check_emass]):
    result = callme\
             (
               structure, escan, vff, direction,
               order = 2,
               nbpoints = 3,
               # will save output to directory "name".
               directory = "Emass_%i_%s" % (i, name),
               # computes at kpoint of deformed structure.
               kpoint = deform_kpoint(kpoint, structure.cell, relaxed.cell),
               # number of states 
               nbstates = nbstates,
               # sets folded method's energy reference
               reference = ref,
             )
    # checks expected are as expected. 
#   assert norm( result[0] - expected ) < 1e-6, "%s\n%s" % (result[0], expected)
    # And print.
    if world.rank == 0:
      print "Ok - %s(%s): %s -> %s: %s" % (name, direction, kpoint, escan.kpoint, result[0])
  exit(0)
  
