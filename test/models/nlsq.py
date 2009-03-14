class Functional:
  def __init__( self, _functional, _structures ):
    from lada import crystal, models, opt

    self.structures = [ crystal.sStructure(s) for s in _structures ]
    self.functional = _functional;

  def __call__( self, _args ):

    from lada import crystal

    self.__argstofunctional__( _args )
    result = float(0)
    for structure in self.structures:
      forces = crystal.sStructure( structure )
#     print structure, forces
      intermed = self.functional( structure, forces ) - structure.energy
      result += intermed * intermed
    return result

  def ngradient( self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=1, stepsize=1e-3 )

  def __argstofunctional__( self, _args ):

    index = 0
    for bond in self.functional.bonds:
      self.functional.bonds[bond.key()].hardsphere = _args[index]; index += 1
      self.functional.bonds[bond.key()].andderwalls = _args[index]; index += 1
    for charge in self.functional.charges:
      self.functional.charges[charge.key()] = float(_args[index]);
      index += 1

  def args( self ):
    from lada import opt
    result = opt.cReals()
    for bond in self.functional.bonds:
      result.append( bond.data().hardsphere )
      result.append( bond.data().vanderwalls )
    result.extend( [ charge.data() for charge in self.functional.charges ] )
    return result

  def __str__( self ):
    string = ""
    for bond in self.functional.bonds:
      string += "%s %12.2f/r^12 - %8.2f/r^6\n"\
                % (bond.key(), bond.data().hardsphere, bond.data().vanderwalls )
    for charge in self.functional.charges:
      string += "%s %8.2f/r\n"\
                 % (charge.key(), charge.data() )
    string += "mesh: %i %i %i\n" % self.functional.mesh
    string += "lj_cutoff: %f \n" % (self.functional.lj_cutoff)
    string += "ewald_cutoff: %f \n" % (self.functional.ewald_cutoff)
    return string

    


