#! /usr/bin/python

import clj_module 

def mean( _list ):

  return reduce( lambda x,y: x + y, _list, 0 ) / len( _list );

def variance( _list ):
  m = mean( _list )
  return reduce( lambda x,y: x + (y - m) *(y -m), _list, 0 ) / len( _list )


class Correlation( clj_module.StructuresFunctional ):

  def __init__( self, _functional, _structures ):
    clj_module.StructuresFunctional.__init__( self, _functional, _structures )

  def __call__( self, args ):

    from lada import crystal
    import math

    
    clj_energies = []
    lda_energies = []
    for structure in self.structures:
      forces = crystal.sStructure( structure )
      clj_energies.append( self.functional( structure, forces ) )
      lda_energies.append( structure.energy )

    clj_mean = mean( clj_energies )
    lda_mean = mean( lda_energies )
    co = 0
    for n, b in enumerate( clj_energies ):
      co += (b - clj_mean) * ( lda_energies[n] - lda_mean )
    co /= math.sqrt( variance( clj_energies ) * variance( lda_energies ) )  * len( clj_energies )
    return co
