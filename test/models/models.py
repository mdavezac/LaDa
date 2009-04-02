#! /usr/bin/python

def conc( _structure, _type ):

  result = 0 
  for atom in _structure.atoms:
    if atom.type == _type: result += 1
  return float( result * 2 ) / float( len( _structure.atoms ) )

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
      func = lambda x: ( str(x[0]), float(x[1]) * 1000.0, 3 )
    else:
      func = lambda x: ( str(x[0]), float(x[1]) * 1000.0, float(x[2]) )
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

def save( _func, _filename ):

  import pickle
  file = open( _filename, "w" )
  pickle.dump( _func, file )
  file.close()

def load( _filename ):
  import pickle
  file = open( _filename, "r" )
  return pickle.load( file )

def reduce_structures( _structures ):

  index = 1
  s1 = _structures[0]
  comp = lambda s0:   abs( s0.energy - s1.energy )\
                    + abs( conc( s0, "Li" ) - conc( s1, "Li" ) ) > 0.0001
  while index != len( _structures ): 
    _structures[index:] = filter( comp, _structures[index:] )
    s1 = _structures[index]
    index += 1
  return _structures  

def increments( _structures, _minimizer, _functional, _n ):

  from lada import opt
  import random

  args = _functional.args()
  args = opt.cReals( [ random.uniform(0,10) for r in args ] )

  for i in range(1, _n+1):
    l = i * len( _structures ) / _n  
    print "size: ", l
    _functional.structures = _structures[:l] 
    _functional.doprint = 0
    _functional.docharges = 1
    _functional.doconstant = 1
    _functional.set_args( args ) 

    print "Iter 0: ", _functional( args )
    _functional.wenergy, _functional.wstress, _functional.wforces = (1,0,0)
    result = _minimizer( _functional, args )
    _functional.set_args( args ) 
    print _functional

  save( _functional, "__negs" )


def main():
  import nlsq
  from lada import models, minimizer, opt, crystal
  import clj_module
  from math import sqrt
  import random
  import statistics
  import ordering

  random.seed()
  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-2, linestep=1e-2, itermax=1000, \
                verbose = 1, strategy="slowest", uncertainties=1, up=1, gradient=1 )

# structures = read_gsgo_history( "LiCsBr_simple" )
  structures = read_gsgo_history( "history.pop_LiCsBr" )
  structures = reduce_structures( structures )
  baseline = opt.ConvexHull()
  baseline.add( ("CsBr", 0, -3.04396490 * 1000 ) )
  baseline.add( ("LiBr", 1, -3.34619915 * 1000 ) )
  ch = opt.ConvexHull(baseline)
  for s in structures:
    ch.add( ( s, conc( s, "Li" ), s.energy ) )
# for s in structures:
#   if s.energy - baseline( conc(s, "Li") ) <= 0e0: s.weight = 10
  structures = filter( lambda s: s.energy - baseline( conc(s, "Li") ) <= 0e0,
                       structures )

# for s in structures:
#   print s, s.energy
#   print

 
# clj = load( "__charges" )
  clj = create_functional()
  func = nlsq.Functional( clj, structures )

# func = load( "__negs2" )
  func.doprint = 0
  func.docharges = 1
  func.doconstant = 1
  increments( structures, minmizer, func, 3)
# args = func.args()
# args = opt.cReals( [ random.uniform(0,10) for r in args ] )
# func.set_args( args ) 
# print func

# print "Iter 0: ", func( args )
# func.wenergy, func.wstress, func.wforces = (1,0,0)
# result = minmizer( func, args )
# func.set_args( args ) 
# print func
# save( func, "__negs" )
  return

  func.doprint = 0
  for setting in [ (1,0,0), (0,1,0), (0,0,1) ]:
    func.wenergy, func.wstress, func.wforces = setting
    print func( args )
  print func

  corr = statistics.Correlation( clj, structures )
  print "correlations: ", corr( args )
  order = ordering.Order( clj, structures )
  print "order: ", order( args )
  print

  

  for structure in structures:
    if structure.energy - baseline( conc(structure, "Li") ) >= 0e0: continue
#   function = clj_module.ScaleFunction( clj, structure) 
#   args = opt.cReals()
#   args.append( sqrt(structure.scale - 0.5) )
#   result = minmizer( function, args )
#   structure.scale = args[0] * args[0]
    forces = crystal.sStructure( structure )
    energy = func.functional( structure, forces ) + func.constant
    print structure.energy - baseline( conc( structure, "Li") ), (structure.energy - energy ) # , forces, "\n"





if __name__ == "__main__":
  main()
