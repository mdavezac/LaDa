obsolete ~~
import clj_module 

class Functional( clj_module.StructuresFunctional ):

  def __init__( self, _functional, _structures ):

    clj_module.StructuresFunctional.__init__( self, _functional, _structures )
    self.wenergy = 1;
    self.wstress = 1;
    self.wforces = 1;
    self.doprint = 0;

  def __call__( self, _args ):

    from lada import crystal

    self.from_array( _args, self.functional )
    result = float(0)
    energy = float(0)
    stress = float(0)
    force = float(0)
    for structure in self.structures:
      forces = crystal.sStructure( structure )
      intermed = self.functional( structure, forces ) + self.constant
      intermed /= float( len( structure.atoms ) )
      intermed -= structure.energy
      if self.wenergy != 0: 
        s = intermed * intermed *  self.wenergy * structure.weight
        result += s
        energy += s
      if self.wstress != 0: 
        s = self.__stress__( forces ) *  self.wstress * structure.weight
        result += s
        stress += s
      if self.wforces != 0: 
        s = self.__forces__( forces ) *  self.wforces * structure.weight
        result += s
        force += s
    if self.doprint != 0:
      print "__call__ ", result, energy, stress, force
    return result

  @staticmethod
  def __stress__( _forces ):
    from lada import atat

    result = 0
    stress = atat.transpose( _forces.cell ) * _forces.cell
    for i in range(0,3):
      for j in range(i,3):
        result += stress[(i,j)]
    return result / 6.0

  @staticmethod
  def __forces__( _forces ):

    result = 0
    for atom in _forces.atoms:
      result += (atom.pos * atom.pos) * 0.3333333333333
    return result / float( len(_forces.atoms ) )



  def __getstate__(self):
    """Return state values to be pickled."""
#   return (  \
#            self.wenergy, self.wstress, self.wforces )
        
    return ( clj_module.StructuresFunctional.__getstate__( self ),\
             self.wenergy, self.wstress, self.wforces )

  def __setstate__(self, state):
    """Restore state from the unpickled state values."""

    ( newstate, self.wenergy, self.wstress, self.wforces ) = state
    clj_module.StructuresFunctional.__setstate__(self, newstate)


