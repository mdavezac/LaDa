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
  from numpy import matrix
  from lada.crystal import fill_structure

  cell = matrix(lattice.cell) * matrix( [ [4,  0, 0.5],
                                          [0,  1,   0],
                                          [0,  0, 0.5] ] )
  structure = fill_structure(cell)
  N = len(structure.atoms)
  for atom in structure.atoms[:N/2]:
    atom.type = structure.lattice.sites[atom.site].type[0]
  for atom in structure.atoms[N/2:]:
    atom.type = structure.lattice.sites[atom.site].type[1]
  structure.atoms[0], structure.atoms[-1] = structure.atoms[1], structure.atoms[0]

  return structure

def deformed_kpoint(ocell, dcell, lcell, kpoint):
  from numpy import dot, matrix
  k = dot( ocell.T,  dot( matrix(lcell).I.T, kpoint) )
  for i in range(3):
    k[i] = k[i] - float( int(k[i]) )
    if k[i] < 0e0: k[i] += 1e0
    if k[i] > 1e0-1e-6: k[i] = 0e0
  return dot(dcell * k)
  


from math import ceil, sqrt
from os.path import join, exists
from numpy import array as np_array, matrix
from boost.mpi import world
from lada.vff import Vff
from lada.escan import Escan, method

input = "input.xml"

# creates lattice
lattice = create_zb_lattice()

# Creates structure
structure = create_structure()

# some kpoints
X = np_array( [0,0,1], dtype="float64" )
G = np_array( [0,0,0], dtype="float64" )
L = np_array( [0.5,0.5,0.5], dtype="float64" )
W = np_array( [1, 0.5,0], dtype="float64" )

# creates vff 
vff = Vff(input, world)
# creates escan
escan = Escan(input, world)

# launch vff and prints output.
relaxed = vff(structure)
print relaxed

# launch pescan for different kpoints.
for kpoint in [G, X, L, W]:
  escan.kpoint = deformed_kpoint(structure.cell, relaxed.cell, lattice.cell, kpoint)
  print kpoint, " -> ", escan.kpoint
  escan.method = method.full_diagonalisation
  escan.vff_inputfile = "atom.config." + str(world.rank)
  vff.print_escan_input(relaxed, escan.vff_inputfile)
  print escan()
  
