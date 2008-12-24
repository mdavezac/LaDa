#! /usr/bin/python
#
# Version: $Id$
#
from LaDa import rMatrix3d, rVector3d, make_rMatrix3d, Lattice, \
                 Structure, Vff, LayeredVff, Escan
import boost.mpi as mpi
from sys import exit

XMLfilename = "input.xml"

lattice = Lattice()
lattice.fromXML( XMLfilename )
vff = LayeredVff()
vff.structure.fromXML( XMLfilename )
structure =  Structure( vff.structure )
vff.fromXML( XMLfilename )
vff.init()
escan = Escan()
escan.set_mpi( boost.mpi.world )
escan.fromXML( XMLfilename )
# print escan.directory



# length = 31
# scales = [ float(x) / length + 5 for x in range( length ) ]
# directions = [ [1,0,0], [0,1,0], [2,1,0], [2,1,1], [3,1,1], [5,1,1] ]
# for direction in directions:
#   vff.direction = rVector3d( direction )
#   vff.structure = Structure( structure )
#   for x in scales:
#     vff.structure.scale = x
#     print x, vff.evaluate() 
#   print "&"
