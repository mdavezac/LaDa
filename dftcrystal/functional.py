__docformat__ = "restructuredtext en"
from .input import AttrBlock
class Functional(AttrBlock):
  """ Wrapper for the CRYSTAL program. """
  def __init__(self):
    """ Creates the crystal wrapper. """
    from .hamiltonian import Dft
    from .basis import BasisSet
    super(Functional, self).__init__()
    self.basis   = BasisSet()
    """ Holds definition of basis functions. """
    self.dft     = Dft()
    """ Holds definition of functional. """
    self.optgeom = AttrBlock(keyword='OPTGEOM')
    """ Holds definition of geometry optimization. """
    self.title   = None
    """ Title of the calculation. 
    
        Overriden by the name of the input structure, if it exists.
    """

  def read_input(self, tree):
    """ Reads file or string with CRYSTAL input. """
    from ..error import IOError
    from .. import CRYSTAL_geom_blocks as starters

    self.title = tree.keys()[0]
    tree = tree[self.title]
    # read optgeom bit.
    found = False
    for starter in starters:
      if starter in tree.keys(): found = True; break
    if found == False:
      raise IOError('Could not find start of input in file.')
    if 'OPTGEOM' in tree[starter].keys():
      self.optgeom.read_input(tree[starter]['OPTGEOM'])

    # read basis set
    if 'BASISSET' in tree.keys(): 
      self.basis.read_input(tree['BASISSET'])

    # read hamiltonian stuff.
    if 'HAMSCF' in tree.keys(): 
      others = AttrBlock()
      others.read_input(tree['HAMSCF'])
      dft = others._crysinput.pop('DFT', None)
      if dft is not None: self.dft.read_input(dft)
      self._crysinput.update(others._crysinput)

  def print_input(self, **kwargs):
    """ Dumps CRYSTAL input to string. """

    # create optgeom part first, since it needs be inserted in the structure
    # bit. Split along lines and remove empty lines at end.
    # if empty, then make it empty.
    optgeom = self.optgeom.print_input(**kwargs).rstrip().split('\n')
    while len(optgeom[-1].rstrip().lstrip()) == 0: optgeom.pop(-1)
    if len(optgeom) == 2: optgeom = []

    result = ''
    if 'structure' in kwargs:
      structure = kwargs['structure']
      # insert name of structure as title.
      if hasattr(structure, 'name'):
        result += str(structure.name).rstrip().lstrip() + '\n'
      # figure out input of structure. remove empty lines.
      lines = structure.print_input(**kwargs).split('\n')
      while len(lines[-1].rstrip().lstrip()) == 0: lines.pop(-1)
      # insert optgeom inside the geometry bit.
      if len(optgeom) > 0: lines.insert(-2, optgeom)
      # turn back into string.
      result += '\n'.join(lines)
    else: # no structures. Not a meaningful input, but whatever.
      result += '\n'.join(optgeom)

    # now add basis
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    result += self.basis.print_input(**kwargs)

    # now add hamiltonian
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    result += self.dft.print_input(**kwargs)

    # add keywords contained directly in functional.
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    a = AttrBlock()
    a._crysinput = self._crysinput
    result += a.print_input(**kwargs)

    # end input and return
    result = result.rstrip()
    if result[-1] != '\n': result += '\n'
    return result + 'END\n'
