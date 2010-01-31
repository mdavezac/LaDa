#! /usr/bin/python
def create_lattice():

  from numpy import array as np_array
  from lada import crystal

  lattice = crystal.Lattice()

  lattice.cell = np_array( [ [    0,  0.5,  0.5 ], \
                             [  0.5,    0,  0.5 ], \
                             [  0.5,  0.5,    0 ] ], dtype="float64" )
  lattice.scale = 4.42

  lattice.sites.append( crystal.Site( np_array([0, 0, 0], dtype="float64"), ["K", "Rb"]) )
  lattice.find_space_group()

  return lattice

# def main():
from lada import ce

lattice = create_lattice()

pair_classes = ce.create_pairs( lattice, 20, 0 )
j = 0
print len(pair_classes)
for pairs in pair_classes:
  j += len(pairs)
  print j
