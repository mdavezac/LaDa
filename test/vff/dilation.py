#! /usr/bin/python
#
# Version: $Id$
#
from LaDa import rMatrix3d, rVector3d, Lattice, \
                 Structure, Vff, LayeredVff
from sys import exit
import boost.mpi as mpi
import pickle 

lattice = Lattice()
lattice.fromXML( "input.xml" )
vff = LayeredVff()
vff.structure.fromXML( "input.xml" )
structure =  Structure( vff.structure )
vff.fromXML( "input.xml" )
vff.init()

length = 31
scales = [ float(x) / length + 5 for x in range( length ) ]
directions = [ [1,0,0] ] #, [0,1,0], [2,1,0], [2,1,1], [3,1,1], [5,1,1] ]
for direction in directions:
  vff.direction = rVector3d( direction )
  vff.structure = Structure( structure )
  if mpi.world.rank == 0:
    vff.evaluate() 
  n = mpi.broadcast( mpi.world, vff.structure, 0 )
  if mpi.world.rank == 1:
    print n
# for x in scales:
#   vff.structure.scale = x
#   print x, vff.evaluate() 
# print "&"
