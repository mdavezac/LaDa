#! /usr/bin/python
import clj_module 

def conc( _structure, _type ):

  result = 0 
  for atom in _structure.atoms:
    if atom.type == _type: result += 1
  return float( result * 2 ) / float( len( _structure.atoms ) )


class Order(clj_module.StructuresFunctional):

  def __init__( self, _functional, _structures, _baseline ):
    from lada import crystal

    clj_module.StructuresFunctional.__init__( self, _functional, _structures )
    self.baseline = _baseline
    efunc = lambda x: x.energy / len(x.atoms) - self.baseline( conc( x, "Li" ) )
    self.structures.sort( lambda x,y:  cmp( efunc(x), efunc(y) ) )
    self.functional = _functional;
    assert len( self.structures ) >= 2, "Number of structures too small."
    self.sigma = self.structures[1].energy - self.structures[0].energy
    for i in range(2, len( self.structures ) ):
      min = self.structures[i].energy - self.structures[i-1].energy
      if self.sigma > min: self.sigma = min
    self.sigma = 10


  def __call__( self, _args ):
    from lada import crystal

    self.from_array( _args, self.functional )
    result = float(0)
    forces = crystal.sStructure( self.structures[0] )
    last =   self.functional( self.structures[0], forces ) \
           / len( self.structures[0].atoms ) \
           - self.baseline( conc( self.structures[0], "Li" ) ) 
    for structure in self.structures[1:]:
      forces = crystal.sStructure( structure )
      intermed =   self.functional( structure, forces ) / len( structure.atoms )\
                 - self.baseline( conc( structure, "Li" ) )
      result += self.__fermi__( float(intermed - last) )
      last = intermed
    print result
    return result

  def __fermi__( self, _x ):
    from math import exp
    a = _x * self.sigma
    if a > 1e3: return 0e0
    if a < -1e3: return 1e0
    c = exp( _x * self.sigma )
    if c > 1e8: return 0e0
    if c < -1e8: return 1e0
    return 1e0 / ( c + 1e0 )




