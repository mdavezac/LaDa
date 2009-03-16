#! /usr/bin/python

class Function:
  def __init__( self, _clj, _structure ):

    from lada import crystal
    self.clj = _clj
    self.structure = crystal.sStructure( _structure )

  def __call__( self, _args ):
    from lada import crystal, models
    models.fold_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    return self.clj( self.structure, forces )

  def ngradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    models.fold_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    self.clj( self.structure, forces )
    self.clj.gradient( self.structure, forces )
    models.unfold_structure( forces, _gradients )


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
  models.unfold_structure( structure, args )

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
# models.fold_structure( args, structure )
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
