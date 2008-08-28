#! /usr/bin/python
from readpi import importPiStructure, SFtoCE_Structure
from LaDa import rMatrix3d, rVector3d, make_rMatrix3d, Lattice, \
                 tetragonalCS
from sys import exit
from numpy import array
from EquivStructure import makeRotationGroup, importLDA,\
                           importStructures


lattice = Lattice()
lattice.fromXML( "input.xml" )

# h = lattice.syms()
# # print "nb syms: ", len( lattice.syms() )
# # for x in lattice.syms():
# #   if x[1] == rVector3d( [0,0,0] ):
# #     print "sym:\n", x[0]

# G = [ x[0] for x in lattice.syms() if x[0].det() > 0 ]

# string = "array( [ [%2i,%2i,%i], [%2i,%2i,%2i], [%2i,%2i,%2i] ] ),"
# for x in G:
#   print  string % ( x[0,0], x[0,1], x[0,2],
#                     x[1,0], x[1,1], x[1,2],
#                     x[2,0], x[2,1], x[2,2] )

cs = tetragonalCS()
cs.fromXML( "input.xml" )
lda = importLDA( filename="originalLDAs.dat" )
for key in lda.keys():
  esStruct = importStructures( key, lda[key], "."  )
  structure = SFtoCE_Structure( esStruct )
  cs.define( structure )
  cs.vars()[:] = [ 2*x-1 for x in esStruct.entries ]
  print "%40s  %6.2f + %6.2f" % \
        ( esStruct.name,
          lda[key] - cs.evaluate(), cs.evaluate() )
  
# filename = "PIs2thru16"
# file = open(filename, 'r')
# oldcell= array( [ [0,0,0], [0,0,0], [0,0,0] ] ) 
# while( True ):

#   try:
#     esStruct = importPiStructure( file )
#   except IOError: 
#     print "End of %s\n" % ( file )
#     exit();
#   else:

#     structure = Structure()
#     if abs( array( oldcell - esStruct.period ).sum() ) > 0.001:
#       structure = SFtoCE_Structure( esStruct )
#       cs.define( structure )
#       oldcell = esStruct.period
#     cs.vars()[:] = [ 2*x-1 for x in esStruct.entries ]
#     print "%8.2f" % ( cs.evaluate() )
#     cs.vars()[:] = [ -x for x in cs.vars() ]
#     print "%8.2f" % ( cs.evaluate() )


#file.close()

