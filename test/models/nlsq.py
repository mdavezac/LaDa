class Functional:
  def __init__( self, _functional, _structures ):
    from lada import crystal, models, opt

    self.structures = [ crystal.sStructure(s) for s in _structures ]
    self.functional = _functional;

  def __call__( self, _args ):

    self.__argstofunctional__( _args )
    result = float(0)
    for structure in self.structures:
      forces = crystal.sStructure( structure )
      intermed = _functional( structure, force ) - structure.energy
      result += intermed * intermed
    return result;

  def ngradient( self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=1, stepsize=1e-3 )

  def __argstofunctional__( self, _args ):

    index = 0
    for bond in self.functional.bonds:
      index += 2
#     bonds[bond.key()].hardsphere = _args[index]; index += 1
#     bonds[bond.key()].vandderwalls = _args[index]; index += 1
    for charge in self.functional.charges:
      self.functional.charges[charge.key()] = float(_args[index]);
      index += 1

  def args( self ):
    from lada import opt
    result = opt.cReals()
    for bond in self.functional.bonds:
      result.append( bond.data().hardsphere )
      result.append( bond.data().vandderwalls )
    result.extend( [ charge.data() for charge in self.functional.charges ] )
    return result

  def __str__( self ):
    string = ""
    for bond in self.functional.bonds:
      string += "%s %12.2f/r^12 - %8.2f/r^6\n"\
                % (bond.key(), bond.data().hardsphere, bond.data().vandderwalls )
    for charge in self.functional.charges:
      string += "%s %8.2f/r\n"\
                 % (charge.key(), charge.data() )
    return string

    


