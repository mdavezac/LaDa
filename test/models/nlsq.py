class Functional:
  def __init__( self, _functional, _structures ):
    from lada import crystal, models, opt

    self.structures = [ crystal.sStructure(s) for s in _structures ]
    self.functional = _functional;
    self.wenerg = 1;
    self.wstress = 1;
    self.wforces = 1;
    self.charges = 0;

  def __call__( self, _args ):

    from lada import crystal

    self.from_array( _args, self.functional )
    result = float(0)
    energy = float(0)
    stress = float(0)
    force = float(0)
    for structure in self.structures:
      forces = crystal.sStructure( structure )
      intermed = self.functional( structure, forces ) - structure.energy
      if self.wstress != 0: 
        s = intermed * intermed *  self.wenergy
        result += s
        energy += s
      if self.wstress != 0: 
        s = self.__stress__( forces ) *  self.wstress
        result += s
        stress += s
      if self.wforces != 0: 
        s = self.__forces__( forces ) *  self.wforces
        result += s
        force += s
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

  def gradient( self, _args, _gradients ):
    from lada import minimizer

    print "starting grad"
    minimizer.interpolated_gradient( self, _args, _gradients,
                                     tolerance=1e-18, n=3, stepsize=1e-3 )
    print "__gradient__", _gradients 
    return _gradients


  def __str__( self ):
    string = "Nlsq fitting functional\n"
    for bond in self.functional.bonds:
      string += "  %s %12.2f/r^12 - %8.2f/r^6\n"\
                % (bond.key(), bond.data().hardsphere, bond.data().vanderwalls )
    for charge in self.functional.charges:
      string += "  %s %8.2f/r\n"\
                 % (charge.key(), charge.data() )
    string += "  mesh: %i %i %i\n" % self.functional.mesh
    string += "  lj_cutoff: %f \n" % (self.functional.lj_cutoff)
    string += "  ewald_cutoff: %f \n" % (self.functional.ewald_cutoff)
    return string

  def from_array( self, _args, _functional ):

    index = 0
    length = len( _args )
    for bond in _functional.bonds:
      _functional.bonds[bond.key()].hardsphere = abs(_args[index]);  index += 1
      _functional.bonds[bond.key()].vanderwalls = abs(_args[index]); index += 1
    sum = 0
    if self.charges == 0: return
    for charge in _functional.charges:
      sign = 1
      if _functional.charges[charge.key()] > 0:
        _functional.charges[charge.key()] = abs(_args[index]);
      else:
        _functional.charges[charge.key()] = -abs(_args[index]);

  def from_functional( self, _functional, _args ):

    index = 0
    length = len( _args )
    for bond in _functional.bonds:
      _args[index] = _functional.bonds[bond.key()].hardsphere;  index += 1
      _args[index] = _functional.bonds[bond.key()].vanderwalls; index += 1
    if self.charges == 0: return
    for charge in _functional.charges:
      _args[index] = abs( _functional.charges[charge.key()] )
      break

  def args(self):
    from lada import opt

    if self.charges != 0: self.charges = 1
    result = opt.cReals( [ 0 for i in range( 0, len(self.functional.bonds) * 2 + self.charges ) ] )
    self.from_functional( self.functional, result )
    return result

  def set_args(self, _args):
    from lada import opt

    self.from_array( _args, self.functional )

