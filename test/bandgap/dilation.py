#! /usr/bin/python
#
# Version: $Id$
#

def nbstates( _escan, _structure ):
  import LaDa
  if _escan.parameters.method == LaDa.full_diagonalization: 
    _escan.parameters.nbstates = LaDa.nb_valence_states( _structure ) + 2
    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      _escan.parameters.nbstates /= 2

def print_bands( _escan, _structure, _first, _last, _result, _nbpoints = 10 ):
  """ Adds eigenvalues to result from linear interpolation
      between kpoint _first to kpoint _last (excluded) using nbpoints.
  """
  
  import LaDa

  # Creates set of points.
  points = [ float(u) * ( _last - _first ) / float(_nbpoints) + _first for u in range(_nbpoints) ]
  for k in points:
    _escan.parameters.kpoint = k
    nbstates( _escan, _structure )
    _escan.run()
    result = [ LaDa.norm2(_escan.parameters.kpoint) ]

    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      for a in _escan.eigenvalues:
        result.append( a ) 
        result.append( a ) 
    else:
      for n in range( len(_escan.eigenvalues)  ):
        result.append( _escan.eigenvalues[ n ] ) 
    _result.append( result )

def main():
  import LaDa
  from LaDa import rMatrix3d, rVector3d, Lattice, \
                   Structure, Vff, LayeredVff, Escan
  import boost.mpi as mpi
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

  vff.evaluate()
  escan.vff_inputfile = "atom_input." + str( mpi.world.rank )
  mpi.broadcast( mpi.world, vff.structure, 0 )
  vff.print_escan_input( escan.vff_inputfile )
  escan.scale = vff.structure

  result = []
  print_bands( escan, vff.structure, rVector3d( [0,0,0] ), rVector3d( [2.0,0,0] ), result, 10 )

  file = open( "result", "w" )
  for eigs in result:
    string = ""
    for u in eigs:
      string += " %f " % ( u )
    print >>file, string

if __name__ == "__main__":
  main()
