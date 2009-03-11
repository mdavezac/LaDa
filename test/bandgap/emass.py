#! /usr/bin/python
#
# Version: $Id$
#
def compute_bandgap( _bandgap, _Evbm, _Ecbm ):
  from lada import atat

  oldnbstates = _bandgap.parameters.nbstates
  if atat.norm2(_bandgap.parameters.kpoint) < 1e-6:
    _bandgap.parameters.nbstates /= 2

  result = []
  _bandgap.parameters.Eref =  _Evbm
  result.extend( sorted( get_values( _bandgap ) ) )
  _bandgap.parameters.Eref =  _Ecbm
  result.extend( sorted( get_values( _bandgap ) ) )
  if atat.norm2(_bandgap.parameters.kpoint) < 1e-6:
    result.extend( result )
    result.sort()
    _bandgap.parameters.nbstates = oldnbstates
  return result

def get_values( _escan ):
  _escan.run()
  return [ float(u) for u in _escan.eigenvalues ]

def read_structure( _xmlinput ): 
  from lada import crystal

  lattice = crystal.Lattice()
  lattice.fromXML( _xmlinput )
  structure = crystal.Structure()
  structure.fromXML( _xmlinput )
  return structure


def compute_distorted_kpoint( _original_cell, _distorted_cell, _kpoint ):
  """
      Computes the position of the kpoint in a distorted lattice, from its 
      position in a non-distorted lattice.
  """

  from lada import atat 

  distorted_kcell = atat.inverse( atat.transpose( _distorted_cell ) )
  return distorted_kcell * atat.transpose( _original_cell ) * _kpoint

def get_masses( _xmlinput, _kpoint, _direction, _step, _nbeval ):

  from lada import crystal, vff, atat, escan

  lattice = crystal.Lattice()
  lattice.fromXML( _xmlinput )
  func_vff = vff.Vff()
  func_vff.structure.fromXML( _xmlinput )
  structure =  crystal.Structure( func_vff.structure )

  original_cell = atat.rMatrix3d( structure.cell )
  dirn = atat.rVector3d( _direction )
  dirn = dirn * 1e0 / sqrt( atat.norm2(dirn) )

  func_vff.fromXML( _xmlinput )
  func_vff.init()
  func_escan = escan.Escan()
  func_escan.set_mpi( mpi.world )
  func_escan.fromXML( _xmlinput )

  func_bg = escan.BandGap()
  func_bg.set_mpi( mpi.world )
  func_bg.fromXML( _xmlinput )

  func_vff.evaluate()
  distorted_cell = atat.rMatrix3d( func_vff.structure.cell )
  bandgap.vff_inputfile = "atomic_config." + str( mpi.world.rank )
  mpi.broadcast( mpi.world, func_vff.structure, 0 )
  func_vff.print_escan_input( bandgap.vff_inputfile )
  bandgap.scale = func_vff.structure
  bandgap.parameters.kpoint = compute_distorted_kpoint\
                              (
                                original_cell,
                                distorted_cell,
                                atat.rVector3d( _kpoint )
                              )
  korigin = atat.rVector3d( bandgap.parameters.kpoint )
  gamma = escan.Bands(bandgap.bands)

  bandgap.escan.Eref = bandgap.Eref.cbm;
  func_emass = escan.eMass()
  result = func_emass\
           (\
             bandgap,
             original_cell, 
             vf.structure,
             _kpoint,
             _direction,
             2,
             bandgap.Eref.cbm
           )

  print result



def create_results( _xmlinput, _kpoint, _direction, _step, _nbval, _filename ):

  import pickle
  from lada import crystal, vff, atat, escan
  import sys
  import boost.mpi as mpi
  from math import sqrt

  lattice = crystal.Lattice()
  lattice.fromXML( _xmlinput )
  func_vff = vff.Vff()
  func_vff.structure.fromXML( _xmlinput )
  structure =  crystal.Structure( func_vff.structure )

  original_cell = atat.rMatrix3d( structure.cell )
  dirn = atat.rVector3d( _direction )
  dirn = dirn * 1e0 / sqrt( atat.norm2(dirn) )

  func_vff.fromXML( _xmlinput )
  func_vff.init()
  bandgap = escan.BandGap()
  bandgap.set_mpi( mpi.world )
  bandgap.fromXML( _xmlinput )

  func_vff.evaluate()
  distorted_cell = atat.rMatrix3d( func_vff.structure.cell )
  bandgap.vff_inputfile = "atomic_config." + str( mpi.world.rank )
  mpi.broadcast( mpi.world, func_vff.structure, 0 )
  func_vff.print_escan_input( bandgap.vff_inputfile )
  bandgap.scale = func_vff.structure
  bandgap.parameters.kpoint = compute_distorted_kpoint\
                              (
                                original_cell,
                                distorted_cell,
                                atat.rVector3d( _kpoint )
                              )
  korigin = atat.rVector3d( bandgap.parameters.kpoint )
  bandgap.evaluate( func_vff.structure )
  gamma = escan.Bands(bandgap.bands)
  bandgap.parameters.nbstates = 4
  bandgap.parameters.method = escan.folded
  result = []

  dummy = [-range(1, _nbval)[-u] for u in range( 1, _nbval )];
  dummy.append( 0 )
  dummy.extend( range(1, _nbval) );
  kpoints = [ 
              compute_distorted_kpoint\
              (
                original_cell,
                distorted_cell,
                  atat.rVector3d( _kpoint )
                + atat.rVector3d( _direction ) * float(_step) * float(u)
              )
              for  u in dummy 
            ]

  # -kpoint
  for k in kpoints:
    bandgap.parameters.kpoint = atat.rVector3d( k )
    dummy = [ bandgap.parameters.kpoint * dirn - atat.rVector3d( korigin ) * dirn ]
    dummy.extend( compute_bandgap( bandgap, gamma.vbm, gamma.cbm ) )
    result.append( dummy )

  file = open( _filename, "w" )
  pickle.dump( (gamma, atat.rVector3d(_kpoint), atat.rVector3d(_direction), result), file )
  file.close()

def read_results( _filename):
  import pickle
  file = open( _filename, "r" )
  result = pickle.load( file)
  file.close
  return result

def print_results( filename ):

  from math import sqrt
  from lada import physics, atat

  (gamma, kpoint, direction, results) = read_results( filename )
  print "# kpoint: ", kpoint,  "  -- direction: ", direction
  for r in results:
    string = "%f" % ( r[0] )
    for u in r[1:]:
      string = "%s %f" % ( string, u )
    print string
# for r in results[steps-1:]:
#   string = "%f" % ( sqrt(LaDa.norm2( r[0] )) )
#   for u in r[1:]:
#     string = "%s %f" % ( string, u / LaDa.Hartree("eV")  )
#   print string
  print "&"

def interpolate_bands( _filename, _scale, _order = 2 ):

  from math import sqrt, pow, pi
  from lada import atat, physics, minimizer
  
  (gamma, kpoint, direction, results) = read_results( _filename )
  print "# kpoint: ", kpoint,  "  -- direction: ", direction
  start = 0
  end = len( results ) - start
  middle = len( results ) / 2 
  for band in range( 1, len( results[0] ) ):
    matrixA = [ \
                [ \
                  pow( r[0], i )\
                  for i in range(0, _order+1)  \
                ]\
                for r in results[start:end] \
              ]
    vectorB = [ r[band] for r in results[start:end] ]
    vectorX = [ 0 for r in range(0, _order+1) ]
    for i in range( 0, len(results) ):
      v = float(abs(i - middle))
      if v == 0: continue
      matrixA[i] = [ u / pow( v, 4 ) for u in matrixA[i] ]
      vectorB[i] /= pow( v, 4) 
    ( x, resid, iter ) = minimizer.linear_lsq( A=matrixA, x=vectorX, b=vectorB, \
                                               verbosity=0, tolerance = 1e-18, itermax = 10000 )
#   print "# ", x, resid, iter
#   for a in results:
#     print a[0], sum( [ x[i]*pow(a[0],i) for i in range(0, _order+1) ] )
#   print "&"
    mass =   physics.Hartree("eV") * 2e0 * pi * pi \
           * physics.a0("A") * physics.a0("A") / _scale / _scale / x[2]
    print "# ", mass

def main():
  # from sys import exit


  scale = read_structure( "sigeemass.xml" ).scale 
  get_masses( "sigeemass.xml", [0,0,0], [1,0,0], 0.01, 10 )
# pickle_filename = "_si.0.01"
# create_results( "sigeemass.xml", [0,0,0], [1,0,0], 0.01, 10, "_ge_gamma" )
# create_results( "sigeemass.xml", [0.5,0.5,0.5], [0.5,0.5,0.5], 0.01, 10, "_ge_Ll" )
# create_results( "sigeemass.xml", [0.5,0.5,0.5], [0.5,-0.5,0], 0.01, 10, "_ge_Lt" )
# create_results( "sigeemass.xml", [1,0,0], [1,0,0], 0.01, 10, "_ge_Xl" )
# print_results( "_ge_gamma" )
# print_results( "_ge_gamma" )
# print_results( "_ge_gamma" )
# print_results( "_ge_Xl" )
# interpolate_bands( "_ge_gamma", scale, 2 )
# print_results( "_si_large_mesh" )
# print_results( pickle_filename )

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
