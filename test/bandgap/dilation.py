#! /usr/bin/python
#
# Version: $Id$
#
from LaDa import rMatrix3d, rVector3d, make_rMatrix3d, Lattice, \
                 Structure, Vff, LayeredVff, Escan
import boost.mpi as mpi
import pickle 
# from sys import exit

XMLfilename = "../input.xml"

lattice = Lattice()
lattice.fromXML( XMLfilename )
vff = LayeredVff()
vff.structure.fromXML( XMLfilename )
structure =  Structure( vff.structure )
vff.fromXML( XMLfilename )
vff.init()
escan = Escan()
escan.set_mpi( mpi.world )
escan.fromXML( XMLfilename )

pickle.dump( vff.structure, open( "pickle", "w" ) )

vff.evaluate()
print vff.structure
# escan.vff_inputfile = "atom_input." + str( mpi.world.rank )
# mpi.broadcast( mpi.world, vff.structure, 0 )
# vff.print_escan_input( escan.vff_inputfile )
# escan.run()
# print escan.eigenvalues



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
