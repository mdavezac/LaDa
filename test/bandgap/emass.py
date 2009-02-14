#! /usr/bin/python
#
# Version: $Id$
#

def get_values( _escan ):
  _escan.run()
  return [ float(u) for u in _escan.eigenvalues ]

def read_structure( _xmlinput ): 
  import LaDa

  lattice = LaDa.Lattice()
  lattice.fromXML( _xmlinput )
  structure = LaDa.Structure()
  structure.fromXML( _xmlinput )
  return structure

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
  dummy = list( result[:len(result)-1] )
  dummy.reverse()
  result.extend( dummy )
# for i in range( 1, _nbval ):
#   bandgap.parameters.kpoint = LaDa.rVector3d( _direction ) * float(_step) * float(i)
#   dummy = [ LaDa.rVector3d( bandgap.parameters.kpoint ) ]
#   bandgap.parameters.Eref =  gamma.vbm
#   dummy.extend( sorted( get_values( bandgap ) ) )
#   bandgap.parameters.Eref =  gamma.cbm
#   dummy.extend( sorted( get_values( bandgap ) ) )
#   result.append( dummy ) 

  file = open( _filename, "w" )
  pickle.dump( (gamma, result), file )
  file.close()

def read_results( _filename):
  import pickle
  file = open( _filename, "r" )
  result = pickle.load( file)
  file.close
  return result

def print_results( filename ):

  from math import sqrt
  import LaDa

  (gamma, results) = read_results( filename )
  steps = len( results ) / 2 + 1 
  for r in results[0:steps-1]:
    string = "%f" % ( -sqrt(LaDa.norm2( r[0] )) )
    for u in r[1:]:
      string = "%s %f" % ( string, u / LaDa.Hartree("eV")  )
    print string
  for r in results[steps-1:]:
    string = "%f" % ( sqrt(LaDa.norm2( r[0] )) )
    for u in r[1:]:
      string = "%s %f" % ( string, u / LaDa.Hartree("eV")  )
    print string
  print "&"

def interpolate_bands( _filename, _scale, _order = 2 ):

  from math import sqrt, pow, pi
  import LaDa
  
  print _scale
  kscale = 2e0 * pi / _scale
  (gamma, results) = read_results( _filename )
  steps = len( results ) / 2 + 1
# for band in range( 1, len( results[0] ) ):
  for band in range( 1, len( results[0] ) ):
    matrixA = [ \
                [ \
                  pow( -sqrt( LaDa.norm2( r[0] ) ) * kscale, i )\
                  for i in range(0, _order+1)  \
                ]\
                for r in results[:steps-1] \
              ]
    matrixA.extend\
    ( \
      [ \
        [ \
          pow( sqrt( LaDa.norm2( r[0] ) ) * kscale , i ) \
          for i in range(0, _order+1)  \
        ]\
        for r in results[steps-1:] \
      ] \
    )
    vectorB = [ r[band] / LaDa.Hartree("eV") for r in results ]
    vectorX = [ 0 for r in range(0, _order+1) ]
    ( x, resid, iter ) = LaDa.linear_lsq( A=matrixA, x=vectorX, b=vectorB, \
                                          verbosity=0, tolerance = 1e-18, itermax = 10000 )
#   print "# ", x, resid, iter
#   for a in range( -steps, steps + 1):
#     u = float(a) / 100.0 
#     print u, sum( [ x[i]*pow(u,i) for i in range(0, _order+1) ] )
#   print "&"
    print "# ", 0.5/x[2]

def main():
  # from sys import exit


  scale = read_structure( "sigeemass.xml" ).scale 
  pickle_filename = "_si.0.01"
# create_results( "sigeemass.xml", [1,0,0], 0.01, 20, pickle_filename )
# print_results( pickle_filename )
# print_results( "_si_large_mesh" )
# print_results( pickle_filename )
  interpolate_bands( pickle_filename, scale, 3 )

#   for i,r in enumerate( matrixA ):
#     print r[1], vectorB[i]
#   print "&"

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
