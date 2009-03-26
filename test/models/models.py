#! /usr/bin/python

def create_functional():

  from math import pow
  from lada import models

  clj = models.Clj()
  clj.mesh = (5, 5, 5)
  clj.lj_cutoff = 10.5
  clj.ewald_cutoff = 25
  clj.charges["Li"] = 1.0
  clj.charges["Cs"] = 1.0
  clj.charges["Br"] = -1.0
  hs = { "Li":1.34, "Cs":2.25,"Br":1.14}
  vdw = {"Li":2.2,"Cs":3,"Br":1.9}
  for a in ["Li", "Cs", "Br" ]:
    for b in ["Li", "Cs", "Br" ]:
      type = models.bond_type( a, b )
      hs_ = float( hs[a] ) + float( hs[b]  )
      vdw_ = float(vdw[a]) + float(vdw[b] )
      clj.bonds[type] = models.LJBond( pow(hs_, 12.0), pow(vdw_, 6.0) )

  return clj

def read_gsgo_history( _filename ):

  def read_structure( _lines ):
  
    from lada import crystal, atat, physics

    while len( _lines ):
      line = _lines.pop(0)
      splitted = line.split()
      if len( splitted ) == 0: continue
      if splitted[0] == "Structure": break
    if not len( _lines ): return (0, 0)

    structure = crystal.sStructure()
    structure.scale = 5

    splitted = _lines.pop(0).split()
    func = 0
    if len( splitted ) < 3:
      func = lambda x: ( str(x[0]), float(x[1]), 3 )
    else:
      func = lambda x: ( str(x[0]), float(x[1]), float(x[2]) )
    (structure.name, structure.energy, structure.scale) = func(splitted)

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
  from lada import models, minimizer, opt, crystal
  import clj_module
  from math import sqrt
  import random

  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-2, linestep=1e-2, itermax=100, \
                verbose = 1, strategy="fast", uncertainties=1, up=1, gradient=1 )

  structures = read_gsgo_history( "LiCsBr_simple" )
  for s in structures:
    print s, s.energy
    print

  epinput = "licsf.input"
  clj = create_functional()
  func = nlsq.Functional( clj, structures )
  func.wenergy = 0 
  func.wstress = 1 
  func.wforces = 1 
  args = func.args()
  args = opt.cReals( [ 1 for r in args ] )
# random.seed()
# args = opt.cReals( [ r * ( 1 + random.uniform(-0.02, 0.02 ) ) for r in args ] )
# func.set_args( args ) 
  print func

  print "Iter 0: ", func( args )
  result = minmizer( func, args )
  func.set_args( args ) 
  print func

  for structure in structures:
  # function = clj_module.ScaleFunction( clj, structure) 
  # args = opt.cReals()
  # args.append( sqrt(structure.scale) )
  # result = minmizer( function, args )
  # structure.scale = args[0] * args[0]
    forces = crystal.sStructure( structure )
    structure.energy = clj( structure, forces )
    print structure.energy, structure.scale, "\n", forces.cell, "\n"





if __name__ == "__main__":
  main()
