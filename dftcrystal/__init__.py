__docformat__ = "restructuredtext en"

from geometry import AffineTransform, AtomicSymmetries, DisplaceAtoms, ExtPrnt,\
                     Ghosts, InsertAtom, MakeSAED, Marker, ModifySymmetry,     \
                     PrnSymDir, RemoveAtom, Slabcut, Slabinfo, Stop, Crystal,\
                     SymmDir, Slab
registered = { 'atomrot':    AffineTransform,
               'atomicsymm': AtomicSymmetries,
               'atomdisp':   DisplaceAtoms,
               'extprnt':    ExtPrnt,
               'ghosts':     Ghosts,
               'atominse':   InsertAtom,
               'makesaed':   MakeSAED,
               'marker':     Marker,
               'modifysymm': ModifySymmetry,
               'prnsymdir':  PrnSymDir,
               'atomremo':   RemoveAtom,
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
