#! /usr/bin/python
#
# Version: $Id$
#

def get_values( _escan ):
  _escan.run()
  return [ float(u) for u in _escan.eigenvalues ]

def create_results( _xmlinput, _direction, _step, _nbval, _filename ):

  import LaDa
  import pickle
  import sys
  import boost.mpi as mpi

  lattice = LaDa.Lattice()
  lattice.fromXML( _xmlinput )
  vff = LaDa.Vff()
  vff.structure.fromXML( _xmlinput )
  structure =  LaDa.Structure( vff.structure )


  vff.fromXML( _xmlinput )
  vff.init()
  bandgap = LaDa.BandGap()
  bandgap.set_mpi( mpi.world )
  bandgap.fromXML( _xmlinput )

  vff.evaluate()
  bandgap.vff_inputfile = "atomic_config." + str( mpi.world.rank )
  mpi.broadcast( mpi.world, vff.structure, 0 )
  vff.print_escan_input( bandgap.vff_inputfile )
  bandgap.scale = vff.structure
  bandgap.evaluate( vff.structure )
  gamma = LaDa.Bands(bandgap.bands)
  bandgap.parameters.nbstates = 4
  bandgap.parameters.method = LaDa.folded
  result = []

  # -kpoint
  for i in [range(1, _nbval)[-u] for u in range( 1, _nbval )]:
    bandgap.parameters.kpoint = LaDa.rVector3d( _direction ) * float(-_step) * float(i)
    dummy = [ LaDa.rVector3d( bandgap.parameters.kpoint ) ]
    bandgap.parameters.Eref =  gamma.vbm
    dummy.extend( sorted( get_values( bandgap ) ) )
    bandgap.parameters.Eref =  gamma.cbm
    dummy.extend( sorted( get_values( bandgap ) ) )
    result.append( dummy )

  # gamma
  bandgap.parameters.nbstates = 2
  bandgap.parameters.kpoint = LaDa.rVector3d( [0,0,0] )
  bandgap.parameters.Eref =  gamma.vbm
  dummy = get_values( bandgap )
  bandgap.parameters.Eref =  gamma.cbm
  dummy.extend( get_values( bandgap ) )
  dummy.extend( dummy )
  dummy.sort()
  dummy.insert( 0, LaDa.rVector3d( bandgap.parameters.kpoint ) )
  result.append( dummy )
  bandgap.parameters.nbstates = 4

  # +kpoint
  for i in range( 1, _nbval ):
    bandgap.parameters.kpoint = LaDa.rVector3d( _direction ) * float(_step) * float(i)
    dummy = [ LaDa.rVector3d( bandgap.parameters.kpoint ) ]
    bandgap.parameters.Eref =  gamma.vbm
    dummy.extend( sorted( get_values( bandgap ) ) )
    bandgap.parameters.Eref =  gamma.cbm
    dummy.extend( sorted( get_values( bandgap ) ) )
    result.append( dummy )

  file = open( _filename, "w" )
  pickle.dump( (gamma, result), file )
  file.close()

def read_results( _filename):
  import pickle
  file = open( _filename, "r" )
  result = pickle.load( file)
  file.close
  return result

def main():
  from math import sqrt
  import LaDa
  # from sys import exit



# create_results( "input.xml", [1,0,0], 0.01, 4, "pickle.bulk" )
  (gamma, results) = read_results( "pickle.bulk" )
# print gamma
  for r in results[0:3]:
    string = "%f" % ( -sqrt(LaDa.norm2( r[0] )) )
    for u in r[1:]:
      string = "%s %f" % ( string, u )
    print string
  for r in results[3:]:
    string = "%f" % ( sqrt(LaDa.norm2( r[0] )) )
    for u in r[1:]:
      string = "%s %f" % ( string, u )
    print string
  print "&"

  for band in range( 3, 4 ): # len( results[0] ) ):
    matrixA = [ [ 1, -sqrt( LaDa.norm2( r[0] ) ), LaDa.norm2( r[0] ) ] for r in results[:3] ]
    matrixA.extend( [ [ 1, sqrt( LaDa.norm2( r[0] ) ), LaDa.norm2( r[0] ) ] for r in results[3:] ] )
    vectorB = [ r[band] for r in results ]
    vectorX = [ 0 for r in range(3) ]
#   print matrixA
#   print vectorB
    ( x, resid, iter ) = LaDa.linear_lsq( A=matrixA, x=vectorX, b=vectorB, \
                                          verbosity=0, tolerance = 1e-18 )
    print resid, iter
    print LaDa.mul_mat_vec( matrixA, x )
    print x
    for a in range( 9 ):
      u = float(a-4) / 100.0 
      print u, x[0] + u * x[1] + u * u * x[2]
    print "&"

#   print A
#   print b
#   print x
#   print 
# print x
# print resid
# print b
# print sum( [ abs(u) for u in LaDa.mul_mat_vec( A, x ) ] )

if __name__ == "__main__":
  main()
