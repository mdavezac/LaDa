from lada.pcm import Clj

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
  structures = structures[:5]

# for s in structures:
#   print s, s.energy
#   print

 
# clj = load( "__charges" )
  clj = create_functional()
  func = nlsq.Functional( clj, structures )
# save( func, "__negs" )

# func = load( "__8" )
  func.structures = structures
  func.doprint = 0
  func.docharges = 1
  func.doconstant = 1
# func = increments( structures, minmizer, func, 3)
# save( _functional, "__negs2" )
# f = load( "__negs" )
# print f
  args = func.args()
# args = opt.cReals( [ random.uniform(0,10) for r in args ] )
  func.set_args( args ) 
# print func

  print "Iter 0: ", func( args )
  func.wenergy, func.wstress, func.wforces = (1,0,0)
# result = minmizer( func, args )
  func.set_args( args ) 
  print func
# save( func, "__8" )
# return

  func.doprint = 0
  for setting in [ (1,0,0), (0,1,0), (0,0,1) ]:
    func.wenergy, func.wstress, func.wforces = setting
    print func( args )
  print func

  corr = statistics.Correlation( clj, structures )
  print "correlations: ", corr( args )
  order = ordering.Order( clj, structures, baseline )
  print "order: ", order( args )
  print
  result = minmizer( order, args )
  order.set_args( args ) 
  print "order: ", order( args )

  

  func.wenergy, func.wstress, func.wforces = (1,0,0)
  for structure in structures:
#   if structure.energy - baseline( conc(structure, "Li") ) >= 0e0: continue
#   function = clj_module.ScaleFunction( clj, structure) 
#   args = opt.cReals()
#   args.append( sqrt(structure.scale - 0.5) )
#   result = minmizer( function, args )
#   structure.scale = args[0] * args[0]
    forces = crystal.sStructure( structure )
    energy =  ( func.functional( structure, forces ) + func.constant ) / len( structure.atoms )
    print structure.energy - baseline( conc( structure, "Li") ),\
         energy - baseline( conc( structure, "Li") )





if __name__ == "__main__":
  main()
