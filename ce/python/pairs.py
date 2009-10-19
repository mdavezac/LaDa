#! /usr/bin/python
def create_lattice():

  from lada import crystal, atat

  lattice = crystal.Lattice()

  lattice.cell = atat.rMatrix3d( [ [    0,  0.5,  0.5 ], \
                                   [  0.5,    0,  0.5 ], \
                                   [  0.5,  0.5,    0 ] ] )
  lattice.scale = 4.42

  lattice.sites.append( crystal.Site( (0, 0, 0) ) )
  lattice.sites[0].type = crystal.StringVector( [ "K", "Rb" ] );
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
