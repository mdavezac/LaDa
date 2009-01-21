#! /usr/bin/python
#
# Version: $Id$
#

def epi_c( _vff ):
  from LaDa import norm2
  from math import sqrt
  u = _vff.direction / sqrt( norm2( _vff.direction ) ) 
  result = u * _vff.structure.cell 
  return result

def outer_prod( _a, _b ):
  from LaDa import rMatrix3d
  return rMatrix3d( [ [ a * b for b in _b ] for a in _a ] )

def delta_dir( _vff, _cell, _dir ):
  from LaDa import norm2
  from math import sqrt
  u = _dir / sqrt( norm2( _dir ) ) 
  result = u * ( _vff.structure.cell - _cell )
  return result



def main():
  from LaDa import rMatrix3d, rVector3d, Lattice, \
                   Structure, Vff, LayeredVff, norm2
  from sys import exit
  import boost.mpi as mpi
  from math import sqrt

  lattice = Lattice()
  lattice.fromXML( "../quaternary.xml" )
  vff = LayeredVff()
  vff.structure.fromXML( "../quaternary.xml" )
  structure =  Structure( vff.structure )
  vff.fromXML( "../quaternary.xml" )
  vff.init()

  # length = 31
  # scales = [ float(x) / length + 5 for x in range( length ) ]
  # directions = [ [1,1,1] ]
  # vff.direction = rVector3d( [0,0,1] )
  # for dir in directions:
  #   vff.structure = structure
  #   vff.direction = rVector3d( dir )
  #   for x in scales :
  #     vff.structure.scale = x
  #     e = vff.evaluate()
  #     print vff.structure.cell * vff.structure.scale
  #     print x, e
  #   print "&"

  length = 30
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  directions = [( rVector3d([1,0,0]), rVector3d([0,1,0]) ), \
                ( rVector3d([2,1,1]), rVector3d([0.5,-0.5,-0.5]) ), \
                ( rVector3d([1,0,1]), rVector3d([1,1,-1]) ), \
                ( rVector3d([1,1,1]), rVector3d([0.5,-0.5,0]) ) ] 
  for dir in directions:
    vff.structure = Structure( structure )
    vff.direction = dir[0]
    print "# ", vff.direction
    for x in scales :
      vff.structure.scale = x
      e = vff.evaluate()
      c = delta_dir( vff, structure.cell, dir[0] )
      a = delta_dir( vff, structure.cell, dir[1] )
      print "%f %s" % (x, e )
    print "&"

if __name__ == "__main__":
  main()
