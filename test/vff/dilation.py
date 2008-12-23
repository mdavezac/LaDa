#! /usr/bin/python
#
# Version: $Id$
#
from LaDa import rMatrix3d, rVector3d, make_rMatrix3d, Lattice, \
                 Structure, Vff, LayeredVff
from sys import exit


lattice = Lattice()
lattice.fromXML( "input.xml" )
vff = LayeredVff()
vff.structure.fromXML( "input.xml" )
structure = vff.structure
vff.fromXML( "input.xml" )
vff.init()

length = 31
scales = [ float(x) / length + 5 for x in range( length ) ]
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 
print "&"

vff.direction = rVector3d( [ 1, 0, 0 ] )
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 
print "&"

vff.direction = rVector3d( [ 0, 1, 0 ] )
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 
print "&"

vff.direction = rVector3d( [ 0, 0, 1 ] )
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 
print "&"

vff.direction = rVector3d( [ 1, 0, 1 ] )
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 
print "&"

vff.direction = rVector3d( [ 2, 0, 1 ] )
vff.structure = structure
for x in scales:
  vff.structure.scale = x
  print x, vff.evaluate() 



