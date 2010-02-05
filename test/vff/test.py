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

from sys import exit
from os.path import join, exists
from numpy import dot
from boost.mpi import world
from lada.vff import Vff

input = "input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates structure and test result
structure, result = create_structure()

# creates vff 
vff = Vff(input, world)

# launch vff and prints output.
relaxed, stress = vff(structure)
if world.rank == 0:
  diff = 0e0
  for a, b in zip(relaxed.atoms, result.atoms):
    assert a.type == b.type
    assert a.site == b.site
    diff += dot(a.pos-b.pos, a.pos-b.pos)

  assert diff / float(len(structure.atoms)) < 1e-8, diff 
  assert abs(relaxed.energy - result.energy) < 1e-8, abs(relaxed.energy - result.energy) 


