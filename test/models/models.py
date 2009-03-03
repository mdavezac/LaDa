#! /usr/bin/python

class Function:
  def __init__( _self, _clj, _structure ):
    _self.clj = _clj
    _self.structure = _structure

  def __call__( _self, _args ):
    from lada import crystal, models
    models.fold_structure( _args, _self.structure )
    forces = crystal.sStructure( _self.structure )
    return _self.clj( _self.structure, forces )

  def gradient( _self, _args, _gradients ):
    from lada import crystal, models, minimizer

    models.fold_structure( _args, _self.structure )
    forces = crystal.sStructure( _self.structure )
    _self.clj( _self.structure, forces )
    _self.clj.gradient( _self.structure, forces )
    models.unfold_structure( forces, _gradients )


  def ngradient( _self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients, n=1, stepsize=1e-3 )

class Parabola:

  def __call__( _self, _args ):
    return _args[0]

  def gradient( _self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients, n=0, stepsize=1e-4, tolerance=1e-18 )

def main():
  from lada import crystal, models, atat, minimizer, opt

  s=str("")
  clj = models.Clj();
  species = []
  epinput = "simple.ep.input"
  poscar = "simple.POSCAR_0"
  models.read_epinput(clj, species, epinput )

  structure = crystal.sStructure();
  crystal.read_poscar( structure, poscar, species )
  crystal.to_fractional( structure );
  structure.atoms[0].freeze = 7;
  forces = crystal.sStructure( structure );
  function = Function( clj, structure) 
  args = opt.cReals()
  models.unfold_structure( structure, args )

  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-0, verbose = 0 )
  result = minmizer( function, args )
# minmizer.set( type="gsl_sd", convergence=1e-12, linestep=1e-0, verbose = 0 )
# result = minmizer( function, args )
# models.fold_structure( args, structure )
  b = clj(structure, forces )
  print "m: ", b, structure, forces
  gradients = opt.cReals( [ 0 for r in args ] )
  ngradients = opt.cReals( [ 0 for r in args ] )
  function.ngradient( args, ngradients )
  function.gradient( args, gradients )
  for i,g in enumerate(ngradients):
    print "%2i %18.9e %18.9e %18.9e" % (i, g, gradients[i], g - gradients[i] )
 
if __name__ == "__main__":
  main()
