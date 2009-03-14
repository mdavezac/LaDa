#! /usr/bin/python

class Function:
  def __init__( self, _clj, _structure ):
    self.clj = _clj
    self.structure = _structure

  def __call__( self, _args ):
    from lada import crystal, models
    models.fold_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    return self.clj( self.structure, forces )

  def gradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    models.fold_structure( _args, self.structure )
    forces = crystal.sStructure( self.structure )
    self.clj( self.structure, forces )
    self.clj.gradient( self.structure, forces )
    models.unfold_structure( forces, _gradients )


  def ngradient( self, _args, _gradients ):
    from lada import crystal, models, minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=1, stepsize=1e-3 )

class Parabola:

  def __call__( self, _args ):
    return _args[0]

  def gradient( self, _args, _gradients ):
    from lada import minimizer

    minimizer.interpolated_gradient( self, _args, _gradients, n=0, stepsize=1e-4, tolerance=1e-18 )

def read_functional( _filename ):
  from lada import crystal, models, atat, minimizer, opt

  clj = models.Clj();
  species = []
  models.read_epinput(clj, species, _filename )

  return clj, species

def read_gsgo_history( _filename ):

  def read_structure( _lines ):
  
    from lada import crystal, atat, physics

    while len( _lines ):
      line = _lines.pop(0)
      if line.split()[0] == "Structure": break
    if not len( _lines ): return (0, 0)

    structure = crystal.sStructure()
    structure.scale = 5

    func = lambda x: (str(x[0]), float(x[1]))
    (structure.name, structure.energy) = func(_lines.pop(0).split())

#   print map( lambda x: float(x), _lines.pop(0).split() )

    structure.cell = atat.rMatrix3d( [ map( lambda x: float(x), _lines.pop(0).split() )
                                       for u in range(0, 3) ] )
    nbatoms = int( _lines.pop(0) )
    _lines.pop(0); _lines.pop(0); _lines.pop(0)
    atomtype = ("Li", "Cs", "Cl" )

    func = lambda x: ( physics.Symbol(int(x[0])),\
                       atat.rVector3d( [ float(x[i]) for i in range(1,4) ] ) )
    for n in range(0, nbatoms):
      atom = crystal.StrAtom()
      ( atom.type, atom.pos ) = func( _lines.pop(0).split() )
      structure.atoms.append( atom )
     
    return (1,structure)


  file=open(_filename,'r')
  lines = file.readlines()
  file.close()
  neof = 1
  structures = []
  while neof:
    (neof, structure) = read_structure( lines )
    if not neof: break
    structures.append( structure )
  return structures


def main():
  import nlsq
  from lada import models, minimizer, opt

  structures = read_gsgo_history( "LiCsBr_simple" )
# for s in structures:
#   print s
 #  print

  epinput = "licsf.input"
# clj, species = read_functional( epinput )
  clj = models.Clj()
  clj.mesh = (3, 3, 3)
  clj.lj_cutoff = 2.5
  clj.ewald_cutoff = 25
  clj.charges["Li"] = 1
  clj.charges["Cs"] = 1
  clj.charges["Br"] = -1
  for bond in ["Li Li", "Li Cs", "Li Br", "Cs Br", "Cs Cs", "Br Br"]:
    clj.bonds[bond] = models.LJBond( 1, 1 )
# nlsq_func = nlsq.Functional( clj, structures )
# print nlsq_func
# args =  nlsq_func.args()

  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-0, verbose = 0 )
  for structure in structures:
    print structure
    function = Function( clj, structure) 
    args = opt.cReals()
    models.unfold_structure( structure, args )
    result = minmizer( function, args )
#   b = clj(structure, forces )
#   print "m: ", b, structure, forces

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
  args = opt.cReals()
  models.unfold_structure( structure, args )

  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-0, verbose = 1, strategy="slowest" )
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
  main2()
