""" Python wrapper for the CRYSTAL dft code. """
__docformat__ = "restructuredtext en"
__all__ = [ 'Extract', 'AffineTransform', 'DisplaceAtoms', 'InsertAtoms',
            'Marker', 'ModifySymmetry', 'RemoveAtoms', 'Slabcut', 'read',
            'Slabinfo','Crystal', 'Molecule', 'Slab', 'Functional', 'Shell'  ]

from .basis import Shell
from .functional import Functional
from .extract import Extract
from .geometry import AffineTransform, DisplaceAtoms, InsertAtoms, Marker,      \
                      ModifySymmetry, RemoveAtoms, Slabcut, Slabinfo, Slab
from .crystal import Crystal
from .molecule import Molecule
registered = { 'atomrot':    AffineTransform,
               'atomdisp':   DisplaceAtoms,
               'atominse':   InsertAtoms,
               'marker':     Marker,
               'modisymm':   ModifySymmetry,
               'atomremo':   RemoveAtoms,
               'slabcut':    Slabcut,
               'slab':       Slab,
               'slabinfo':   Slabinfo,
               'crystal':    Crystal,
               'molecule':   Molecule }
""" Keywords used in creating the geometry. """

def read(path):
  """ Reads functional and structure from stream """
  from .parse import parse
  b = parse(path)
  key = b.keys()[0]

  func = Functional()
  crys = Crystal()
  func.read_input(b)
  crys.read_input(b[key]['CRYSTAL'])
  crys.name = key
  return crys, func
