#! /usr/bin/python
#
# Version: $Id$
#

def get_values( _escan ):
  _escan.run()
  return [ float(u) for u in _escan.eigenvalues ]

def create_results( _bandgap, _structure, _direction, _step, _nbval, _filename ):

  import LaDa
  import pickle
  import sys

  _bandgap.evaluate( _structure )
  gamma = LaDa.Bands(_bandgap.bands)
  _bandgap.parameters.Eref =  ( gamma.cbm + gamma.vbm ) * 0.5
  _bandgap.parameters.nbstates = 4
  _bandgap.parameters.method = LaDa.folded
  result = []
  for i in range( 1, _nbval ):
    _bandgap.parameters.kpoint = -_direction * float(_step) * float(i)
    result.append( ( _bandgap.parameters.kpoint, list( get_values( _bandgap ) ) ) )
  for i in range( 0, _nbval ):
    _bandgap.parameters.kpoint = _direction * float(_step) * float(i)
    result.append( ( _bandgap.parameters.kpoint, list( get_values( _bandgap ) ) ) )

  file = open( "pickle", "w" )
  pickle.dump( result, file )
  file.close()

def read_results( _filename):
  import pickle
  file = open( _filename, "r" )
  result = pickle.load( file)
  file.close
  return result

def main():
  import LaDa
  from LaDa import rMatrix3d, rVector3d, Lattice, \
                   Structure, LayeredVff, BandGap, Bands
  import boost.mpi as mpi
  # from sys import exit


  XMLfilename = "input.xml"

  lattice = Lattice()
  lattice.fromXML( XMLfilename )
  vff = LayeredVff()
  vff.structure.fromXML( XMLfilename )
  structure =  Structure( vff.structure )


  vff.fromXML( XMLfilename )
  vff.init()
  bandgap = BandGap()
  bandgap.set_mpi( mpi.world )
  bandgap.fromXML( XMLfilename )

  vff.evaluate()
  bandgap.vff_inputfile = "atomic_config." + str( mpi.world.rank )
  mpi.broadcast( mpi.world, vff.structure, 0 )
  vff.print_escan_input( bandgap.vff_inputfile )
  bandgap.scale = vff.structure

# create_results( bandgap, vff.structure, rVector3d( [1,0,0] ), 0.01, 1, "pickle" )
  results = read_results( "pickle" )
  print results[0][0], results[0][1]

if __name__ == "__main__":
  main()
