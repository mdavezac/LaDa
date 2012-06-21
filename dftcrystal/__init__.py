__docformat__ = "restructuredtext en"

from geometry import AffineTransform, AtomicSymmetries, DisplaceAtoms, ExtPrnt,\
                     Ghosts, InsertAtoms, MakeSAED, Marker, ModifySymmetry,    \
                     PrnSymDir, RemoveAtoms, Slabcut, Slabinfo, Stop, Crystal, \
                     SymmDir, Slab
registered = { 'atomrot':    AffineTransform,
               'atomsymm':   AtomicSymmetries,
               'atomdisp':   DisplaceAtoms,
               'extprnt':    ExtPrnt,
               'ghosts':     Ghosts,
               'atominse':   InsertAtoms,
               'makesaed':   MakeSAED,
               'marker':     Marker,
               'modisymm':   ModifySymmetry,
               'prnsymdir':  PrnSymDir,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'stop':       Stop,
               'structure':  Crystal,
               'symmdir':    SymmDir, }

# class Structure(object):
#   """ Functional structure class. 
#   
#       Forms of the basis of the functional approach CRYSTAL_ uses to define
#       structures.

#       .. _CRYSTAL: http://www.crystal.unito.it/
#   """
#   def __init__(self, spacegroup='

def call_crystal(input, program=None):
  """ Performs call to crystal in current directory """
  pass
