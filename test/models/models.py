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
    _self.clj( _self.structure, forces )
    models.unfold_structure( forces, _gradients )

  def ngradient( _self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients, n = 6, stepsize=1e-6 )

class Parabola:

  def __call__( _self, _args ):
    return _args[0]

  def gradient( _self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( _self, _args, _gradients )

def main():
  from lada import crystal, models

  s=str("")
  clj = models.Clj();
  species = []
  models.read_epinput(clj, species, "ep.input" )

  structure = crystal.sStructure();
  crystal.read_poscar( structure, "POSCAR_0", species )
  forces = crystal.sStructure( structure );
  crystal.to_fractional( structure );
  print clj( structure, forces ), structure, forces
# print structure, forces

  function = Function( clj, structure) 
  args = []
  models.unfold_structure( structure, args )
  print function( args )
  gradients = []
  ngradients = [ 0 for r in args ]
  function.ngradient( args, ngradients )
  function.gradient( args, gradients )
  for i,g in enumerate(ngradients):
    print i, g, gradients[i]

  print 
  print 

  p = Parabola()
  print p([0]), p([1]), p([2])
  ngrad = [0]
  p.gradient([1], ngrad)
  print ngrad[0]

if __name__ == "__main__":
  main()
