#! /usr/bin/python
#
# Version: $Id$
#


def dir_check( _vff, _direction ):
  
  from LaDa import norm2
  from math import sqrt

  u = _direction / sqrt( norm2( _direction ) ) 
  return u * _vff.structure.cell  * _vff.structure.scale


def diffdir_check( _vff, _cell, _direction ):
  
  from LaDa import norm2
  from math import sqrt

  u = _direction / sqrt( norm2( _direction ) ) 
  return u * ( _cell - _vff.structure.cell ) * _vff.structure.scale

def get_asize( _vff, _cell, _direction ):
  
  from LaDa import norm2, inv_rMatrix3d
  from math import sqrt

  u = inv_rMatrix3d( _cell ) * _direction;
  u = _vff.structure.cell * _vff.structure.scale * u;
  return sqrt( norm2(u) / norm2( _direction ) )


def main():
  from LaDa import rMatrix3d, rVector3d, Lattice, \
                   Structure, Vff, LayeredVff, norm2, inv_rMatrix3d
  from sys import exit
  import boost.mpi as mpi
  from math import sqrt

  lattice = Lattice()
  lattice.fromXML( "input.xml" )
  vff = LayeredVff()
  vff.structure.fromXML( "input.xml" )
  structure =  Structure( vff.structure )
  vff.fromXML( "input.xml" )
  vff.init()

  length = 50
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  directions = [ ( rVector3d([1,0,0]), rVector3d([0,1,0]) ), \
                 ( rVector3d([1,1,0]), rVector3d([0, 0, 1]) ), \
                 ( rVector3d([1,1,1]), rVector3d([-1,1,0]) ) ]
  for dir in directions:
    vff.direction = dir[0] 
    print "# ", vff.direction
    for x in scales:
      vff.structure = Structure( structure )
      vff.structure.scale = x
      e = vff.evaluate()
      c = dir_check( vff, dir[0] )
      a = diffdir_check( vff, structure.cell, dir[1] )
      get_asize( vff, structure.cell, dir[1] )
      print get_asize( vff, structure.cell, dir[1] ) / 0.529177, \
            get_asize( vff, structure.cell, dir[0] ) / 0.529177
    print "&"

if __name__ == "__main__":
  main()
