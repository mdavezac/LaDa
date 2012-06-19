from .input import Block, Keyword

class Structure(Block):
  """ CRYSTAL-wise structure, e.g. functional approach. """
  keyword = 'crystal'
  """ CRYSTAL keyword. """
  def __init__(self, spacegroup, *args, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    super(Structure, self).__init__()
    self.spacegroup = spacegroup
    """ Space-group. """
    self.params = args
    """ Lattice parameters. """
    self.shift = kwargs.pop('shift', None)
    """ Shift of the origin of the cell,
    
         - if None, then no shift (equivalent to IFSO == 0)
         - if 0 or second, second setting given international tables
           (equivalent to IFSO == 0).
         - if 1 or first, second setting given international tables
           (equivalent to IFSO == 1).
    """
    self.operators = kwargs.pop('operators', [])
    """ Functional operations on the structure. """
    self.atoms = []
    """ List of atoms. """
    self.shift = None
    for key, value in kwargs.iteritems(): setattr(self, key, value)

  @property
  def shift(self):
    """ Shift of the origin. """
    return self._shift
  @shift.setter
  def shift(self, value):
    from ..error import ValueError
    if value is None: self._shift = None
    elif hasattr(value, '__iter__'): 
      self._shift = array(value, dtype='float64')
    elif value == 0 or value == 'second': self._shift = 0
    elif value == 1 or value == 'first': self._shift = 1
    raise ValueError('Unknown input value for shift.')

  def add_atom(self, *args, **kwargs):
    """ Adds an atom to the structure. """
    from ..crystal import Atom
    atom = Atom(*args, **kwargs)
    self._atoms.append(atom)
    return self

  def __repr__(self, indent=''):
    """ Dumps representation to string. """
    from inspect import getargspec
    # prints construction part.
    length = len(indent + '{0.__class__.__name__}('.format(self))
    result = indent + '{0.__class__.__name__}({0.spacegroup!r}, '.format(self)
    for o in self.params: result += '{0!r}, '.format(o)
    for key, value in self.__dict__.iteritems():
      if key == '_shift': key = 'shift'
      result += '\\\n{length}{0}={1!r}, '.format(key, value, length=length)
    result = result[:-2] + ')'
    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent + '    ', dummy)
    # prints inner keywords/blocks
    if len(self) > 0:
      indent += '    '
      for item in self: 
        hasindent = getargspec(item.__repr__).keywords
        hasindent = False if hasindent is None else 'indent' in hasindent
        result += '\n{0}{1!r}'.format(indent, item) if not hasindent           \
                  else '\n{0}{1}'.format(indent, item.__repr__(indent)) 
    if len(self) > 0:
      result = result[:result.rfind(')')+1]
    return result

  @property 
  def raw(self):
    """ Raw input for structure. """
    from ..physics.spacegroup import space_group as sg
    from ..error import ValueError
    from .. import periodic_table as pt

    if self.shift is None:                hasshift = 0
    elif hasattr(self.shift, '__iter__'): hasshift = 2
    else:                                 hasshift = 1

    if isinstance(self.spacegroup): 
      if self.spacegroup not in sg:
        raise ValueError('Unknown space group {0.spacegroup}'.format(self))
      spacegroup = self.spacegroup
    else: 
      spacegroup = None
      for key, value in sg.iteritems():
        if value.index == self.spacegroup: 
          spacegroup = value.index
          break
      if spacegroup is None:
        raise ValueError('Unknown space group {0.spacegroup}'.format(self))
    
    # control flags + spacegroup 
    result = '1 {0} {1}\n{2}\n'.format( getattr(self, 'ifhr', 0),              \
                                      hasshift, spacegroup )
    # shift if needed
    if hasshift == 2:
      result += '{self.shift[0]} {self.shift[1]} {self.shift[2]}\n'            \
                .format(self)

    # lattice parameters.
    for o in self.params: result += '{0} '.format(o)

    # number of atoms + atoms
    result += '\n{0}\n'.format(len(self))
    for atom in self:
      n = getattr(pt, atom.type, None)
      if n is not None: n = n.atomic_number
      else:
        try: n = int(atom.type)
        except: 
          raise ValueError( 'Could not make sense of atomic type {0.type}.'    \
                            .format(atom) )
      result += '{0} {1.pos[0]} {1.pos[1]} {1.pos[2]}\n'.format(n, atom)
    return result

class GeomKeyword(Keyword):
  """ Adds breaksymm to :py:class:`~lada.dftcrystal.input.Keyword`. """
  def __init__(self, keyword=None, raw=None, *kwargs):
    """ Creates a geometry keyword. 
    
        :param str keyword: 
          keyword identifying the block.
        :param bool keepsym: 
          Whether to keep symmetries or not. Defaults to True.
        :param bool breaksym:
          Whether to break symmetries or not. Defaults to False.
          Only one of breaksymm needs be specified. If both are, then they
          should be consistent.
    """
    super(GeomKeyword, self).__init__(keyword=keyword, raw=raw)
    if 'keepsym' in kwargs and 'breaksym' in kwargs:
      if kwargs['keepsym'] == kwargs['breaksym']:
        raise ValueError('keepsym and breaksym both specified and equal.')
    self.breaksym = kwargs.get('breaksym', not kwargs.get('keepsym', False))
    """ Whether or not to break symmetries.
    
        Defaults to False. If symmetries are not broken, then all equivalent
        atoms are removed.
    """
    kwargs.pop('breaksym', None)
    kwargs.pop('keepsym', None)
  @property
  def keepsym(self):
    """ Not an alias for breaksym. """
    return not self.breaksym
  @keepsym.setter
  def keepsym(self, value):
    self.breaksym = not value
    
class RemoveAtom(GeomKeyword):
  """ Remove atoms from structure. """
  keyword = 'atomremo'
  """ CRYSTAL keyword. """
  def __init__(self, *args, **kwargs):
    """ Creates atom removal operator.
    
        :param args:
          Atoms to remove are set by specifying them as arguments to the
          iterator.
        :param bool breaksym:
          Keyword argument specifying whether to keep or remove symmetries.
          If symmetries are kept, then all equivalent atoms are removed from
          the structure.
        :param bool keepsym:
          Opposite of breaksym. If both are specified, then they should make
          sense...
    """ 
    from copy import copy
    super(RemoveAtom, self).__init__(**kwargs)
    self.labels = copy(args)

  @property
  def raw(self, **kwargs):
    """ Returns input string. """
    return str(len(self.labels)) + '\n'.join(str(u) for u in self.labels) + '\n'

  def __repr__(self):
    result = '{0.__class__.__name__}({1}, '.format(self, ', '.join(self.labels))
    if self.breaksym == True: result += 'keepsym=True, '
    return result[:-2] + ')'
    
class Ghosts(RemoveAtom): 
  """ Sets ghosts atom in structure. """
  keyword = 'ghosts'
  """ CRYSTAL keyword. """
  def __init__(self, *args, **kwargs):
    """ See :py:class:`RemoveAtom` """
    super(Ghosts, self).__init__(*args, **kwargs)
   
class Slabinfo(Keyword):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl):
    from numpy import array
    super(Slabinfo, self).__init__()
    self.hkl = array(hkl, dtype='int32')
    """ Surface normal. """
  @property
  def raw(self):
    """ Input to SLABINFO """
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}'
  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}])'     \
           .format(self)

class Slabcut(Slabinfo):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl, isup, depth):
    super(Slabinfo, self).__init__(hkl)
    self.isup = isup
    """ Subsurface atomic label. """
    self.isup = isup
    """ Subsurface atomic label. """
  @property
  def raw(self):
    """ Input to SLABINFO """
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}\n'
  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}])'     \
           .format(self)

class DisplaceAtoms(GeomKeyword):
  """ Displaces atoms. """
  keyword = 'atomdisp'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    """ Creates a displacement field. """
    super(DisplaceAtoms, self).__init__(**kwargs)
    self.displacements = []
    """ Atomic displacements. """

  def add_atom(self, *args, **kwargs):
    """ Adds a displacement to a given atom.
    
        At present, atom.type should be an index to an atom in the structure.
    """
    from ..crystal import Atom
    self.displacements.append(Atom(*args, **kwargs))

  def __repr__(self, indent= ''):
    """ Dumps representation to string. """
    result = super(DisplaceAtoms, self).__repr__()
    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent + '    ', dummy)
    return result

  @property
  def raw(self):
    """ Raw input to CRYSTAL. """
    result = '{0}\n'.format(len(self.displacements))
    for atom in self.displacements:
      result += '{0.type} {0.pos[0]} {0.pos[1]} {0.pos[2]}\n'.format(atom)
    return result
   
class InsertAtom(GeomKeyword):
  """ insert atom into structure. """
  keyword = 'atominse'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    super(DisplaceAtoms, self).__init__(**kwargs)
    self.atoms = []
    """ Atoms to insert. """

  def add_atom(self, *args, **kwargs):
    """ Adds a displacement to a given atom.
    
        At present, atom.type should be an index to an atom in the structure.
    """
    from ..crystal import Atom
    self.atoms.append(Atom(*args, **kwargs))

  def __repr__(self, indent= ''):
    """ Dumps representation to string. """
    result = super(DisplaceAtoms, self).__repr__()
    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent + '    ', dummy)
    return result

  @property
  def raw(self):
    """ Raw input to CRYSTAL. """
    from .. import periodic_table as pt
    # number of atoms + atoms
    result = '{0}\n'.format(len(self.atoms))
    for atom in self:
      n = getattr(pt, atom.type, None)
      if n is not None: n = n.atomic_number
      else:
        try: n = int(atom.type)
        except: 
          raise ValueError( 'Could not make sense of atomic type {0.type}.'    \
                            .format(atom) )
      result += '{0} {1.pos[0]} {1.pos[1]} {1.pos[2]}\n'.format(n, atom)
    return result

class ModifySymmetry(Keyword):
  """ Modify symmetries. """ 
  keyword = 'modisymm'
  """ CRYSTAL keyword. """
  def __init__(self, *args):
    """ Creates symmetry modifier. 

        Each argument is a sequence of atomic label. Each argument will be
        assigned a different flag. This input is somewhat more condensed than
        the original CRYSTAL input.
    """ 
    super(ModifySymmetry, self).__init__()
    self.atoms = []
    """ Atoms for which to modify symmetries. """
    for o in args: self.atoms.append(o if hasattr(o, '__iter__') else [o])

  @property
  def raw(self):
    """ Raw CRYSTAL input. """
    result = '{0}\n'.format(sum(len(o) for o in self.atoms))
    for labels in self.atoms:
      result += ' '.join(str(u) for u in labels) + '\n'
    return result

class Stop(Keyword):
  """ Stop execution. """
  keyword = 'stop'
  """ CRYSTAL keyword. """

class Marker(Keyword):
  """ Places a marker in the execution list. """
  def __init__(self, name):
    """ Creates a marker. """
    self.keyword = name
    """ Name of marker. """
  @property
  def name(self):
    """ Name of marker. """
    return self.keyword
  @name.setter
  def name(self, value):
    self.keyword = value
  def __repr__(self):
    return '{0.__class__.__name__}({0.name)'.format(self)
  def print_input(**kwargs): return None

class AtomicSymmetries(Keyword):
  """ Prints atomic symmetries. """
  keyword = 'atomsymm'
  """ CRYSTAL keyword. """

class MakeSAED(Keyword):
  """ Prints symmetry allowed distortions. """
  keyword = 'makesaed'
  """ CRYSTAL keyword. """

class PrnSymDir(Keyword):
  """ Prints symmetry allowed displacements. """
  keyword = 'prnsymdir'
  """ CRYSTAL keyword. """

class SymmDir(Keyword):
  """ Prints symmetry allowed geometry optimization directions. """
  keyword = 'symmdir'
  """ CRYSTAL keyword. """

class ExtPrnt(Keyword):
  """ Prints external file format. """
  keyword = 'extprnt'
  """ CRYSTAL keyword. """
