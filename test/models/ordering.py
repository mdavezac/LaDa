#! /usr/bin/python
import clj_module 

class Order(clj_module.StructuresFunctional):

  def __init__( self, _functional, _structures ):
    from lada import crystal

    clj_module.StructuresFunctional.__init__( self, _functional, _structures )
    self.structures.sort( lambda x,y: x.energy < y.energy )
    self.functional = _functional;
    assert len( self.structures ) >= 2, "Number of structures too small."
    self.sigma = self.structures[1].energy - self.structures[0].energy
    for i in range(2, len( self.structures ) ):
      min = self.structures[i].energy - self.structures[i-1].energy
      if self.sigma > min: self.sigma = min
    self.sigma /= 2e0


  def __call__( self, _args ):
    from lada import crystal

    self.from_array( _args, self.functional )
    result = float(0)
    forces = crystal.sStructure( self.structures[0] )
    last = self.functional( self.structures[0], forces )
    for structure in self.structures[1:]:
      forces = crystal.sStructure( structure )
      intermed = self.functional( structure, forces ) / len( structure.atoms )
      result += self.__fermi__( float(intermed - last) )
      last = intermed
    return result

  def __fermi__( self, _x ):
    from math import exp
    return 1e0 / ( exp( _x / self.sigma ) + 1e0 )




