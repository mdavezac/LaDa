#! /usr/bin/python
#
# Version: $Id$
#


def dir_check( _vff, _direction ):
  import numpy as np

  u = _direction / np.linalg.basic.norm(_direction)
  return np.dot(_vff.structure.cell, u)  * _vff.structure.scale


def diffdir_check( _vff, _cell, _direction ):
  
  import numpy as np
  from math import sqrt

  u = _direction / np.linalg.basic.norm(_direction)
  return np.dot(_cell - _vff.structure.cell, u) * _vff.structure.scale

def get_asize( _vff, _cell, _direction ):
  
  import numpy as np
  from math import sqrt

  u = np.dot(_cell.I, _direction);
  u = np.dot(_vff.structure.cell, u) * _vff.structure.scale;
  return np.linalg.basic.norm(u) / np.linalg.basic.norm(_direction) 


def main():
  from sys import exit
  from math import sqrt
  import numpy as np
  import boost.mpi as mpi
  from lada import vff, crystal

  lattice = crystal.Lattice()
  lattice.fromXML( "input.xml" )
  functional = vff.LayeredVff()
  functional.structure.fromXML( "input.xml" )
  structure =  crystal.Structure( functional.structure )
  functional.fromXML( "input.xml" )
  functional.init()

  length = 50
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  directions = [ ( np.array([1,0,0],dtype="float64"), np.array([0,1,0],dtype="float64") ), \
                 ( np.array([1,1,0],dtype="float64"), np.array([0, 0, 1],dtype="float64") ), \
                 ( np.array([1,1,1],dtype="float64"), np.array([-1,1,0],dtype="float64") ) ]
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
