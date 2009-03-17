#! /usr/bin/python

def structure_to_array( _structure, _args ):

  from lada import atat

  _args.clear();
  g = atat.transpose( _structure.cell ) * _structure.cell
  for i in range(0, 3):
    for j in range(i, 3 ):
      _args.append( 0.5 * ( g[(i,j)] + g[(j,i)] ) )
  for atom in _structure.atoms:
    _args.append(atom.pos[0]) 
    _args.append(atom.pos[1]) 
    _args.append(atom.pos[2]) 
  print _args

def array_to_structure( _args, _structure ):

  from lada import atat
  from math import sqrt

  assert len( _args ) == 6 + len( _structure.atoms ) * 3, "wrong size."

  index = 6
  for atom in _structure.atoms:
    atom.pos[0] = _args[index]; index += 1
    atom.pos[1] = _args[index]; index += 1
    atom.pos[2] = _args[index]; index += 1
  a = ( atat.rVector3d(), atat.rVector3d(), atat.rVector3d() )
  for j in range( 0, 3 ):
    for i in range( 0, 3 ):
      a[j][i] = _structure.cell[(j,i)]

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
    print self.structure, result
    return result

  def ngradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    array_to_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    self.clj( self.structure, forces )
    self.clj.gradient( self.structure, forces )
    array_to_structure( forces, _gradients )


  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=3, stepsize=1e-2 )

class ScaleFunction:
  def __init__( self, _clj, _structure ):
    from lada import crystal
    self.clj = _clj
    self.structure = crystal.sStructure(_structure)

  def __call__( self, _args ):
    from lada import crystal, models
    self.structure.scale = _args[0]
    forces = crystal.sStructure( self.structure )
    return self.clj( self.structure, forces )

  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=1, stepsize=1e-3 )

class Parabola:

  def __call__( self, _args ):
    return _args[0]

  def gradient( self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=0, stepsize=1e-4, tolerance=1e-18 )

def main2():
  from lada import crystal, models, atat, minimizer, opt

  s=str("")
  epinput = "ep.input"
  poscar = "POSCAR_0"
  clj, species = read_functional( epinput )

  structure = crystal.sStructure();
  crystal.read_poscar( structure, poscar, species )
  crystal.to_fractional( structure );
  structure.atoms[0].freeze = 7;
  forces = crystal.sStructure( structure );
  function = Function( clj, structure) 

  a = clj(structure, forces )
  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-0, verbose = 0, strategy="slowest" )

  scfunc = ScaleFunction( clj, structure )
  scargs = opt.cReals()
  scargs.append( structure.scale )
  args = opt.cReals()
  array_to_structure( structure, args )

  minmizer( scfunc, scargs )
  result = minmizer( function, args )
  minmizer( scfunc, scargs )
  result = minmizer( function, args )
  minmizer( scfunc, scargs )
  result = minmizer( function, args )
  minmizer( scfunc, scargs )
  result = minmizer( function, args )
  minmizer( scfunc, scargs )
# minmizer.set( type="gsl_sd", convergence=1e-12, linestep=1e-0, verbose = 0 )
# result = minmizer( function, args )
# structure_to_array( args, structure )
  b = clj(structure, forces )
  print "m: ", a, b, structure, forces
  gradients = opt.cReals( [ 0 for r in args ] )
  ngradients = opt.cReals( [ 0 for r in args ] )
  function.ngradient( args, ngradients )
  function.gradient( args, gradients )
  for i,g in enumerate(ngradients):
    print "%2i %18.9e %18.9e %18.9e" % (i, g, gradients[i], g - gradients[i] )
 
if __name__ == "__main__":
  main2()
