#! /usr/bin/python

def structure_to_array( _structure, _args ):

  from lada import atat, models
  from math import acos

  models.unfold_structure( _structure, _args )
  return
  _args.clear();
  g = atat.transpose( _structure.cell ) * _structure.cell
  for i in range(0, 3):
    for j in range(i, 3 ):
      _args.append( 0.5 * ( g[(i,j)] + g[(j,i)] ) )
  for atom in _structure.atoms:
    _args.append( atom.pos[0] )
    _args.append( atom.pos[1] )
    _args.append( atom.pos[2] )

def array_to_structure( _args, _structure ):

  from lada import atat, models
  from math import sqrt, cos

  models.fold_structure( _args, _structure )
  return
  assert len( _args ) == 6 + len( _structure.atoms ) * 3, "wrong size."

  index = 6
  for atom in _structure.atoms:
    atom.pos[0] = _args[index]; index += 1
    atom.pos[1] = _args[index]; index += 1
    atom.pos[2] = _args[index]; index += 1
  a = ( atat.rVector3d(), atat.rVector3d(), atat.rVector3d() )
  for j in range( 0, 3 ):
    for i in range( 0, 3 ):
      a[i][j] = _structure.cell[(i,j)]

  # u, n, m vectors, with u in a0 direction, n perp to a0 with a1 component,...
  u = a[0] / sqrt( atat.norm2(a[0]) )
  n = a[1] - ( a[1] * u ) * u
  n *= (1e0/sqrt( atat.norm2(n) ))
  m = a[2] - ( a[2] * u ) * u - ( a[2] * n ) * n
  m *= 1 / sqrt( atat.norm2(m) )
  unm = atat.rMatrix3d( [ [ u[i], n[i], m[i] ] for i in range(0,3) ] )

  # sign for direct basis. we use det AB = det A det B.
  sign = 1.e0
  if atat.det( unm ) < 0: sign = -1.e0

  # transformation matrix
  gtransform = atat.rMatrix3d()
  assert _args[0] > 0.0001, "null or negative length."
  assert _args[3] > 0.0001, "null or negative length."
  assert _args[5] > 0.0001, "null or negative length."
  gtransform[(0,0)] = sqrt(_args[0])
  gtransform[(0,1)] = _args[1] / gtransform[(0,0)]
  gtransform[(0,2)] = _args[2] / gtransform[(0,0)] 
  assert _args[3] > gtransform[(0,1)] * gtransform[(0,1)], "error." 
  gtransform[(1,1)] = sqrt( _args[3] - gtransform[(0,1)] * gtransform[(0,1)] )
  assert abs( gtransform[(1,1)] ), "error."
  gtransform[(1,2)] = ( _args[4] - gtransform[(0,1)] * gtransform[(0,2)] ) / gtransform[(1,1)]
  func = lambda x,y,z: x - y*y - z*z
  assert func(_args[5], gtransform[(1,2)], gtransform[(0,2)] ) > 0e0, "error"
  func = lambda x,y,z: sqrt( x - y*y - z*z  )
  gtransform[(2,2)] = sign * func( _args[5], gtransform[(1,2)], gtransform[(0,2)] )

  # final result.
  _structure.cell = unm * gtransform


class Function:
  def __init__( self, _clj, _structure ):

    from lada import crystal
    self.clj = _clj
    self.structure = crystal.sStructure( _structure )

  def __call__( self, _args ):
    from lada import crystal, models

    array_to_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    result = self.clj( self.structure, forces )
    return result

  def ngradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer, atat

    array_to_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    self.clj.gradient( self.structure, forces )
    structure_to_array( forces, _gradients )
    return _gradients


  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=2, 
                                     tolerance=1e-12, itermax=100, stepsize=1e-5 )
    return _gradients


class ScaleFunction:
  def __init__( self, _clj, _structure ):
    from lada import crystal
    self.clj = _clj
    self.structure = crystal.sStructure(_structure)

  def __call__( self, _args ):
    from lada import crystal, models
    self.structure.scale = _args[0] * _args[0] + 0.5
    forces = crystal.sStructure( self.structure )
    print "__call__", self.structure.scale
    return self.clj( self.structure, forces )

  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=2, stepsize=1e-1, tolerance=1e-12 )
    return _gradients

class CellFunction:

  def __init__( self, _clj, _structure ):
    from lada import crystal
    self.clj = _clj
    self.structure = crystal.sStructure(_structure)

  def __call__( self, _args ):

    from lada import crystal, models

    self.from_array( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    result = self.clj( self.structure, forces )
    return result

  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    self.from_array( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    self.clj.gradient( self.structure, forces )
    self.from_structure( forces, _gradients )
    return _gradients
    
  @staticmethod
  def from_structure( _structure, _args ):

    for i in range(0, 3):
      for j in range(0, 3):
        _args[ 3*i + j ] = _structure.cell[(i,j)]

  @staticmethod
  def from_array( _array, _structure ):

    for i in range(0, 3):
      for j in range(0, 3):
        _structure.cell[(i,j)] = _array[ 3*i+j]

  @staticmethod
  def args( _structure ):
    from lada import opt

    result = opt.cReals( [ 0 for r in range(0,9) ] )
    CellFunction.from_structure( _structure, result )
    return result

class PosFunction( CellFunction ):

  def __init__( self, _clj, _structure ):
    from lada import crystal
    CellFunction.__init__( self, _clj, _structure )

  @staticmethod
  def from_structure( _structure, _args ):

    index = 0
    for at in _structure.atoms:
      for i in range(0, 3):
        _args[index] = at.pos[i]
        index += 1

  @staticmethod
  def from_array( _array, _structure ):

    index = 0
    for at in _structure.atoms:
      for i in range(0, 3):
        at.pos[i] = _array[index] 
        index += 1

  @staticmethod
  def args( _structure ):
    from lada import opt

    result = opt.cReals( [ 0 for r in range(0, len(_structure.atoms)*3) ] )
    PosFunction.from_structure( _structure, result )
    return result

class StructuresFunctional:

  def __init__( self, _functional, _structures ):
    from lada import crystal, models, opt

    self.structures = [ crystal.sStructure(s) for s in _structures ]
    self.functional = _functional;
    self.doprint = 0
    self.doconstant = 0;
    self.constant = 0;
    self.docharges = 0;

  def gradient( self, _args, _gradients ):
    from lada import minimizer

    self.doprint = 0
    minimizer.interpolated_gradient( self, _args, _gradients,
                                     tolerance=1e-18, n=3, stepsize=1e-3 )
    self.doprint = 1
    return _gradients

  def __str__( self ):
    string = "Nlsq fitting functional\n"
    for bond in self.functional.bonds:
      string += "  %s %12.2f/r^12 - %8.2f/r^6\n"\
                % (bond.key(), bond.data().hardsphere, bond.data().vanderwalls )
    for charge in self.functional.charges:
      string += "  %s %8.2f/r\n"\
                 % (charge.key(), charge.data() )
    if self.doconstant: string += "  constant: %f \n" % self.constant
    string += "  mesh: %i %i %i\n" % self.functional.mesh
    string += "  lj_cutoff: %f \n" % (self.functional.lj_cutoff)
    string += "  ewald_cutoff: %f \n" % (self.functional.ewald_cutoff)
    return string

  def from_array( self, _args, _functional ):

    index = 0
    if self.doconstant:
      index = 1
      self.constant = _args[0] 
    else: self.constant = 0

    for bond in _functional.bonds:
      _functional.bonds[bond.key()].hardsphere = abs(_args[index]);  index += 1
      _functional.bonds[bond.key()].vanderwalls = abs(_args[index]); index += 1
    sum = 0
    if self.docharges == 0: return
    for charge in _functional.charges:
      sign = 1
      if _functional.charges[charge.key()] > 0:
        _functional.charges[charge.key()] = abs(_args[index]);
      else:
        _functional.charges[charge.key()] = -abs(_args[index]);

  def from_functional( self, _functional, _args ):

    index = 0
    if self.doconstant:
      index = 1 
      _args[0] = self.constant

    for bond in _functional.bonds:
      _args[index] = _functional.bonds[bond.key()].hardsphere;  index += 1
      _args[index] = _functional.bonds[bond.key()].vanderwalls; index += 1
    if self.docharges == 0: return
    for charge in _functional.charges:
      _args[index] = abs( _functional.charges[charge.key()] )
      break

  def args(self):
    from lada import opt

    if self.docharges != 0: self.docharges = 1
    result = opt.cReals( [ 0 for i in range( 0, len(self.functional.bonds) * 2 + self.docharges ) ] )
    if self.doconstant: result.append( 0 )
    self.from_functional( self.functional, result )
    return result

  def set_args(self, _args):
    from lada import opt

    self.from_array( _args, self.functional )

  def __getstate__(self):
    """Return state values to be pickled."""

    return ( self.structures, self.functional, self.doprint, \
             self.doconstant, self.constant, self.docharges )

  def __setstate__(self, state):
    """Restore state from the unpickled state values."""

    ( self.structures, self.functional, self.doprint, \
      self.doconstant, self.constant, self.docharges ) = state

class Parabola:

  def __call__( self, _args ):
    from math import cos, sin
    result = _args[0] * _args[0] - cos( _args[0] )
    return result

  def gradient( self, _args, _gradients ):
    from lada import minimizer

#   _gradients[0] = 2* _args[0]
    minimizer.interpolated_gradient( self, _args, _gradients, n=2, stepsize=1e-1, tolerance=1e-12 )
    return _gradients;

def read_functional( _filename ):
  from lada import crystal, models, atat, minimizer, opt

  clj = models.Clj();
  species = []
  models.read_epinput(clj, species, _filename )

  return clj, species


def main2():
  from lada import crystal, models, atat, minimizer, opt

  s=str("")
  epinput = "xiuwen.input"
  clj, species = read_functional( epinput )
  clj.mesh = (8, 8, 8)
  clj.lj_cutoff = 55
  clj.ewald_cutoff = 55

  structure = crystal.sStructure();
  crystal.read_poscar( structure, "POSCAR_0", species )
  crystal.to_fractional( structure );
  function = Function( clj, structure )
  
  rocksalt = crystal.sStructure();
  crystal.read_poscar( rocksalt, "POSCAR_0", species )
  crystal.to_fractional( rocksalt );

  multiple = crystal.sStructure();
  crystal.read_poscar( multiple, "POSCAR_1", species )
  crystal.to_fractional( multiple );

  forcesa = crystal.sStructure( rocksalt )
  a = clj.lennard_johnes( rocksalt, forcesa )
  print "rocksalt: ", a
  forcesb = crystal.sStructure( multiple )
  b = clj.lennard_johnes( multiple, forcesb ) / 2e0
  print "multiple: ", b
  print a, b, a-b
  print forcesa, forcesb

# a = clj(structure, forces )
# minmizer = minimizer.Minimizer()
# minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-0, verbose = 0, strategy="slowest" )

# scfunc = ScaleFunction( clj, structure )
# scargs = opt.cReals()
# scargs.append( structure.scale )
  args = opt.cReals()
  structure_to_array( structure, args )

# minmizer( scfunc, scargs )
# result = minmizer( function, args )
# minmizer( scfunc, scargs )
# result = minmizer( function, args )
# minmizer( scfunc, scargs )
# result = minmizer( function, args )
# minmizer( scfunc, scargs )
# result = minmizer( function, args )
# minmizer( scfunc, scargs )
# minmizer.set( type="gsl_sd", convergence=1e-12, linestep=1e-0, verbose = 0 )
# result = minmizer( function, args )
# structure_to_array( args, structure )
# b = clj(structure, forces )
  gradients = opt.cReals( [ 0 for r in args ] )
  ngradients = opt.cReals( [ 0 for r in args ] )
  function.ngradient( args, ngradients )
  function.gradient( args, gradients )
  for i,g in enumerate(ngradients):
    print "%2i %18.9e %18.9e %18.9e" % (i, g, gradients[i], g - gradients[i] )
 
if __name__ == "__main__":
  main2()
