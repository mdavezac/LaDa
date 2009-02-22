#! /usr/bin/python
#
# Version: $Id$
#


def dir_check( _vff, _direction ):
  
  from lada import atat
  from math import sqrt

  u = _direction / sqrt( atat.norm2( _direction ) ) 
  return u * _vff.structure.cell  * _vff.structure.scale


def diffdir_check( _vff, _cell, _direction ):
  
  from lada import atat
  from math import sqrt

  u = _direction / sqrt( atat.norm2( _direction ) ) 
  return u * ( _cell - _vff.structure.cell ) * _vff.structure.scale

def get_asize( _vff, _cell, _direction ):
  
  from lada import atat
  from math import sqrt

  u = atat.inverse( _cell ) * _direction;
  u = _vff.structure.cell * _vff.structure.scale * u;
  return sqrt( atat.norm2(u) / atat.norm2( _direction ) )


def main():
  from lada import atat, vff, crystal
  from sys import exit
  import boost.mpi as mpi
  from math import sqrt

  lattice = crystal.Lattice()
  lattice.fromXML( "input.xml" )
  functional = vff.LayeredVff()
  functional.structure.fromXML( "input.xml" )
  structure =  crystal.Structure( functional.structure )
  print functional.structure.lattice();
  functional.fromXML( "input.xml" )
  return;
  functional.init()

  length = 50
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  directions = [ ( atat.rVector3d([1,0,0]), atat.rVector3d([0,1,0]) ), \
                 ( atat.rVector3d([1,1,0]), atat.rVector3d([0, 0, 1]) ), \
                 ( atat.rVector3d([1,1,1]), atat.rVector3d([-1,1,0]) ) ]
  for dir in directions:
    functional.direction = dir[0] 
    print "# ", functional.direction
    for x in scales:
      functional.structure = crystal.Structure( structure )
      functional.structure.scale = x
      e = functional.evaluate()
      c = dir_check( functional, dir[0] )
      a = diffdir_check( functional, structure.cell, dir[1] )
      get_asize( functional, structure.cell, dir[1] )
      print get_asize( functional, structure.cell, dir[1] ) / 0.529177, \
            get_asize( functional, structure.cell, dir[0] ) / 0.529177
    print "&"

if __name__ == "__main__":
  main()
