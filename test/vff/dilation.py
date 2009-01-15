#! /usr/bin/python
#
# Version: $Id$
#
from LaDa import rMatrix3d, rVector3d, Lattice, \
                 Structure, Vff, LayeredVff
from sys import exit
import boost.mpi as mpi

lattice = Lattice()
lattice.fromXML( "input.xml" )
vff = LayeredVff()
vff.structure.fromXML( "input.xml" )
structure =  Structure( vff.structure )
vff.fromXML( "input.xml" )
vff.init()

length = 31
scales = [ float(x) / length + 5 for x in range( length ) ]
directions = [ [0,0,1] ]
vff.direction = rVector3d( [0,0,1] )
# for dir in directions:
#   vff.structure = structure
#   vff.direction = rVector3d( dir )
#   for x in scales :
vff.structure.scale = 5.48
e = vff.evaluate()
u = vff.structure.cell[(2,2)]
print x, e, u
# print "&"
