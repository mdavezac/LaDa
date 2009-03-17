#! /usr/bin/python

def create_functional():

  from math import pow
  from lada import models

  clj = models.Clj()
  clj.mesh = (3, 3, 3)
  clj.lj_cutoff = 2.5
  clj.ewald_cutoff = 25
  clj.charges["Li"] = 1.5
  clj.charges["Cs"] = 0.5
  clj.charges["Br"] = -1
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

  structures = read_gsgo_history( "LiCsBr_simple" )
# for s in structures:
#   print s
 #  print

  epinput = "licsf.input"
# clj, species = read_functional( epinput )
  clj = create_functional()
  print clj
# clj = models.Clj()
# clj.mesh = (3, 3, 3)
# clj.lj_cutoff = 2.5
# clj.ewald_cutoff = 25
# clj.charges["Li"] = 1
# clj.charges["Cs"] = 1
# clj.charges["Br"] = -1
# for bond in ["Li Li", "Li Cs", "Li Br", "Cs Br", "Cs Cs", "Br Br"]:
#   clj.bonds[bond] = models.LJBond( 0.5, 1 )
# nlsq_func = nlsq.Functional( clj, structures )
# print nlsq_func
# args =  nlsq_func.args()

  function = clj_module.Parabola()
  args = opt.cReals()
  args.append(5)
  minmizer = minimizer.Minimizer()
  minmizer.set( type="minuit2", convergence=1e-12, linestep=1e-2, itermax=50, \
                verbose = 1, strategy="slowest", uncertainties=1, up=1 )
# minmizer( function, args )
# print args, function(args)
# return
  for structure in structures:
    function = clj_module.ScaleFunction( clj, structure) 
    args = opt.cReals()
    args.append( float(structure.scale) )
    result = minmizer( function, args )
    structure.scale = args[0]

    forces = crystal.sStructure( structure )
    print "energy: ", clj( structure, forces )

    forces = crystal.sStructure( structure )
    a = clj( structure, forces )
    print structure.scale, structure, forces
    function = clj_module.Function( clj, structure) 
    args = opt.cReals()
    clj_module.structure_to_array( structure, args )
    result = minmizer( function, args )
    print result, args
    clj_module.array_to_structure( args, structure )
    forces = crystal.sStructure( structure )
    a = clj( structure, forces )
    print structure.scale, structure, forces
    print "f energy: ", a
#   b = clj(structure, forces )
#   print "m: ", b, structure, forces





if __name__ == "__main__":
  main()
