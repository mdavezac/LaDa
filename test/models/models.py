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
  from lada import models

  structures = read_gsgo_history( "history.pop_LiCsBr" )
# for s in structures:
#   print s
 #  print

  epinput = "ep.input"
  clj, species = read_functional( epinput )
  nlsq_func = nlsq.Functional( clj, structures )
  print nlsq_func

  args = nlsq_func.args()
  print args
  args[0] = 3
  print args

  clj.bonds["Cl Na"] = models.LJBond()
  clj.charges["V"] = 6
  nlsq_func = nlsq.Functional( clj, structures )
  print nlsq_func, clj.charges["V"]

  

# for bond in clj.bonds:
#   print   "%s %8.4f/r^12 - %8.4f/r^6"\
#         % (bond.key(), bond.data().hardsphere, bond.data().vandderwalls )
# for charge in clj.charges:
#   print   "%s %8.4f/r"\
#         % (charge.key(), charge.data() )

def main2():
  from lada import crystal, models, atat, minimizer, opt

  s=str("")
  epinput = "simple.ep.input"
  poscar = "simple.POSCAR_0"
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
