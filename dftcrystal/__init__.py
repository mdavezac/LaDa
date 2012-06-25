__docformat__ = "restructuredtext en"

from geometry import AffineTransform, DisplaceAtoms, Ghosts, InsertAtoms,    \
                     Marker, ModifySymmetry, RemoveAtoms, Slabcut, Slabinfo, \
                     Crystal, Slab
from hamiltonian import Dft, Hybrid, Exchange, Correlation, NonLocal
registered = { 'atomrot':    AffineTransform,
               'atomdisp':   DisplaceAtoms,
               'ghosts':     Ghosts,
               'atominse':   InsertAtoms,
               'marker':     Marker,
               'modisymm':   ModifySymmetry,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'structure':  Crystal,
               'dft':        Dft,
               'exchange':   Exchange,
               'correlat':   Correlation,
               'nonlocal':   NonLocal,
               'hybrid':     Hybrid }

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
