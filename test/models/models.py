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
    from lada import crystal, models

    models.fold_structure( _args, _self.structure )
    forces = crystal.sStructure( _self.structure )
    _self.clj.gradient( _self.structure, forces )
    models.unfold_structure( forces, _gradients )

    return

  def ngradient( _self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients, n = 3, stepsize=1e-2 )

class Parabola:

  def __call__( _self, _args ):
    return _args[0]

  def gradient( _self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients, n=0, stepsize=1e-4 )

def main():
  from lada import crystal, models, atat

  s=str("")
  clj = models.Clj();
  species = []
  models.read_epinput(clj, species, "ep.input" )

  structure = crystal.sStructure();
  crystal.read_poscar( structure, "POSCAR_0", species )
  forces = crystal.sStructure( structure );
  crystal.to_fractional( structure );
  function = Function( clj, structure) 
  args = []
  models.unfold_structure( structure, args )
  print function( args )
  gradients = []
  ngradients = [ 0 for r in args ]
  function.ngradient( args, ngradients )
  function.gradient( args, gradients )
  for i,g in enumerate(ngradients):
    print i, 2* (g - gradients[i]) / ( g+ gradients[i])
 
if __name__ == "__main__":
  main()
