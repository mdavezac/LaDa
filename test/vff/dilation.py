#! /usr/bin/python
#
# Version: $Id$
#


def dir_check( _vff, _direction ):
  from numpy import dot
  from numpy.linalg import norm

  u = _direction / norm(_direction)
  return dot(_vff.structure.cell, u)  * _vff.structure.scale


def diffdir_check( _vff, _cell, _direction ):
  from math import sqrt
  from numpy import dot
  from numpy.linalg import norm

  u = _direction / norm(_direction)
  return dot(_cell - _vff.structure.cell, u) * _vff.structure.scale

def get_asize( _vff, _cell, _direction ):
  from math import sqrt
  from numpy import dot, matrix
  from numpy.linalg import norm

  u = dot( matrix(_cell).I, _direction).T;
  u = dot(_vff.structure.cell, u) * _vff.structure.scale;
  return norm(u) / norm(_direction) 


def main():
  from sys import exit
  from math import sqrt
  from numpy import array
  import boost.mpi as mpi
  from lada import vff, crystal

  lattice = crystal.Lattice("ternary.xml")
  functional = vff.LayeredVff("ternary.xml", mpi.world)
  functional.structure.fromXML( "ternary.xml" )
  structure =  crystal.Structure( functional.structure )
  functional.init()

  length = 50
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  directions = [ ( array([1,0,0],dtype="float64"), array([0,1,0],dtype="float64") ), \
                 ( array([1,1,0],dtype="float64"), array([0, 0, 1],dtype="float64") ), \
                 ( array([1,1,1],dtype="float64"), array([-1,1,0],dtype="float64") ) ]
  for dir in directions:
    functional.direction = dir[0] 
    if mpi.world.rank == 0: print "# ", functional.direction
    for x in scales:
      functional.structure = crystal.Structure( structure )
      functional.structure.scale = x
      e = functional.evaluate()
      c = dir_check( functional, dir[0] )
      a = diffdir_check( functional, structure.cell, dir[1] )
      get_asize( functional, structure.cell, dir[1] )
      if mpi.world.rank == 0:
        print get_asize( functional, structure.cell, dir[1] ) / 0.529177, \
              get_asize( functional, structure.cell, dir[0] ) / 0.529177
    if mpi.world.rank == 0: print "&"

if __name__ == "__main__":
  main()
