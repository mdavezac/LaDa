#! /usr/bin/python
#
# Version: $Id$
#

def nbstates( _escan, _structure ):
  """
     Sets the number of states to be computed in _escan
     to the number of valence band + one conduction band 
     ( + two if spin polarized ) in _structure.
  """
  import LaDa
  if _escan.parameters.method == LaDa.full_diagonalization: 
    _escan.parameters.nbstates = LaDa.nb_valence_states( _structure ) + 2
    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      _escan.parameters.nbstates /= 2

def get_asize( _distorted_structure, _undistorted_cell, _direction ):
  """
     Returns the lattice parameter in _direction in a.u.
  """
  
  from LaDa import norm2, inv_rMatrix3d
  from math import sqrt

  u = inv_rMatrix3d( _cell ) * _direction;
  u = _structure.cell * _structure.scale * u;
  return sqrt( norm2(u) / norm2( _direction ) )

def print_bands( _escan, _structure, _ocell, _first, _last, _result, _nbpoints = 10, _offset = 0 ):
  """
     Adds eigenvalues to _result from linear interpolation
     between kpoint _first to kpoint _last (included) using nbpoints.
     The results are appended to _result as a list of lists:
       [
         [ norm of (undistorted) kpoint, eig1, eig2, ..., eign ],
         [ ....  next kpoint ... ],
         ...
       ]
  """
  
  import LaDa
  from math import sqrt

  # Creates set of points.
  points = [ float(u) * ( _last - _first ) / float(_nbpoints) + _first
             for u in range(_nbpoints+1) ]
  for k in points:
    _escan.parameters.kpoint = compute_distorted_kpoint( _ocell, _structure.cell, k )
    nbstates( _escan, _structure )
    _escan.run()
    result = [ sqrt( LaDa.norm2(k) ) + _offset ]

    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      for a in _escan.eigenvalues:
        result.append( a ) 
        result.append( a ) 
    else:
      for n in range( len(_escan.eigenvalues)  ):
        result.append( _escan.eigenvalues[ n ] ) 
    _result.append( result )

def compute_distorted_kpoint( _original_cell, _distorted_cell, _kpoint ):
  """
      Computes the position of the kpoint in a distorted lattice, from its 
      position in a non-distorted lattice.
  """

  import LaDa

  distorted_kcell = LaDa.inv_rMatrix3d( LaDa.trans_rMatrix3d( _distorted_cell ) )
  return distorted_kcell * LaDa.trans_rMatrix3d( _original_cell ) * _kpoint

def save_strings_to_file( _filename, _strings ):
  """
      Saves a list of strings to file.
  """
  file = open( _filename, "w" )
  for string in _strings:
    print >>file, string 
  file.close()

def compute_allbands( _mpi, _vff, _escan, _structure, _kpoint, _scale, _perpdir, _add_states ):
  """ 
     Computes all valence bands + 1 conduction band + _add_states
     for a set kpoint, scale and direction.
  """

    import LaDa

    # Computes epi structure.
    _vff.structure = LaDa.Structure( _structure )
    _vff.structure.scale = _scale
    epi_energy = _vff.evaluate()

    # computes escan parameters.
    _escan.parameters.kpoint \
       = compute_distorted_kpoint( _structure.cell, _vff.structure.cell, _kpoint )
    nbstates( _escan, _vff.structure )
    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      _escan.parameters.nbstates += _add_states / 2
    else:
      _escan.parameters.nbstates += _add_states 
    _escan.vff_inputfile = "atom_input." + str( _mpi.world.rank )
    _mpi.broadcast( _mpi.world, _vff.structure, 0 )
    _vff.print_escan_input( _escan.vff_inputfile )
    _escan.scale = _vff.structure

    # computes eigs
    _escan.run()

    # prints eigs to a string.
    string = "%f " % ( get_asize( _vff.structure, _structure.cell, _perpdir ) / 0.529177 )

    if    _escan.parameters.potential != LaDa.spinorbit \
       or LaDa.norm2( _escan.parameters.kpoint ) < 1e-6:
      for eig in _escan.eigenvalues: string = "%s %f %f " % ( string,  eig, eig  )
    else:
      for eig in _escan.eigenvalues: string = "%s %f " % ( string, eig )

    return (string, epi_energy)
      
def compute_bandgap( _mpi, _vff, _bandgap, _structure, _kpoint, _scale, _perpdir ):
  """ 
     Computes the bandgap for a set kpoint, scale and direction.
  """

    import LaDa

    # Computes epi structure.
    _vff.structure = LaDa.Structure( _structure )
    _vff.structure.scale = _scale
    epi_energy = _vff.evaluate()

    # computes bandgap parameters.
    _bandgap.vff_inputfile = "atom_input." + str( _mpi.world.rank )
    _mpi.broadcast( _mpi.world, _vff.structure, 0 )
    _vff.print_escan_input( _bandgap.vff_inputfile )
    _bandgap.scale = _vff.structure
    _bandgap.parameters.kpoint \
       = compute_distorted_kpoint( _structure.cell, _vff.structure.cell, _kpoint )

    # computes eigs
    _bandgap.evaluate( _vff.structure )
    # prints eigs to string.
    string =    "%f %s %s "\
             % ( get_asize( _vff.structure, _structure.cell, _perpdir ) / 0.529177,\
                 _bandgap.bands.vbm, _bandgap.bands.cbm )

    return (string, epi_energy)

def main():
  import LaDa
  import boost.mpi as mpi
  from math import sqrt
  # from sys import exit

  # input file from which to read vff, escan, and structure.
  XMLfilename = "dilation.xml"

  # Loads all parameters from file.
  lattice = LaDa.Lattice()
  lattice.fromXML( XMLfilename )
  vff = LaDa.LayeredVff()
  vff.structure.fromXML( XMLfilename )
  structure =  LaDa.Structure( vff.structure )
  vff.fromXML( XMLfilename )
  vff.init()
  escan = LaDa.BandGap()
  escan.set_mpi( mpi.world )
  escan.fromXML( XMLfilename )

  # Defines a few standard kpoint.
  X = LaDa.rVector3d( [1,0,0] )
  G = LaDa.rVector3d( [0,0,0] )
  L = LaDa.rVector3d( [0.5,0.5,0.5] )

  # Defines a few parameters.
  vff.direction = LaDa.rVector3d( [1,1,0] )
  perpdir = LaDa.rVector3d( [-1,1,0] )
  length = 5
  scales = [ 5.43 + x*(5.65-5.43)/ length  for x in range( length + 1 ) ]
  kpoints = ( X, G, L )
  add_states = 6

  results=[]

  # This part prints to file the eigenvalues with respect to kpoint and scale.
  for kpoint in kpoints:

    results.append( "%s %s %s %s" % ( "# epi:", vff.direction, "  kpoint: ", kpoint ) )

    # if computing bandgap only, uncomment this section.
    # if computing all bands, comment it out.
    # begin section.
    escan.parameters.method = LaDa.full_diagonalization
    ( string, epi_energy ) = compute_bandgap( mpi, vff, escan, structure, \
                                              kpoint, scales[0], perpdir )
    results.append( string )
    escan.eref = LaDa.Bands( escan.bands )
    escan.parameters.method = LaDa.folded
    # end section.

    for x in scales[1:]:

      # computes all bands.
      # ( string, epi_energy ) = compute_allbands( mpi, vff, escan, structure, \
      #                                            kpoint, x, perpdir, add_states )
      # computes band gap only
      ( string, epi_energy ) = compute_bandgap( mpi, vff, escan, structure, \
                                                kpoint, x, perpdir )

      results.append( string )

    # end of loop over scales
    results.append( "&" )

  # end of loop over kpoints

  save_strings_to_file( "ge_bands", results )

  # This part computes a band structure for a fixed scale.
# compute_allbands( mpi, vff, escan, structure, \
#                   G, scales[3], perpdir, add_states )
# result = []
# print_bands( escan, vff.structure, structure.cell,
#              LaDa.rVector3d( [0,0,2] ), LaDa.rVector3d( [0,0,1] ),
#              result, 2 )
# file = open( "distorted_periodic", "w" )
# for eigs in result:
#   string = ""
#   for u in eigs:
#     string += " %f " % ( u )
#   print >>file, string


# Tells Python to call function main on execution.
if __name__ == "__main__":
  main()
