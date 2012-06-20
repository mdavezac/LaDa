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
    # shift if needed, watch out for weird ass input.
    if hasshift == 2:
      shift = self.shift * 24.0 + 0.1
      shift = [int(u) for u in shift]
      result += '{0[0]} {0[1]} {0[2]}\n'.format(shift)

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

  @raw.setter
  def raw(self, value):
    """ Reads crystal input. """
    from .. import periodic_table as pt
    from numpy import array
    value = value.split('\n')
    data = value.pop(0).split()
    self.spacegroup, self.ifhr, self.ifso = (int(u) for u in data[:3]) 
    data = value.pop(0)
    if self.spacegroup == 0: 
      self.spacegroup = int(data.split())
    else: self.spacegroup = data.rstrip().lstrip()

    if self.ifso > 1: 
      self.shift = array(value.pop(0).split(), dtype='float64') / 24.0

    self.params = value.pop(0).split()
    for line in value:
      line = line.split()
      type = int(line[0])
      if type < 100: 
        for key, value in pt.__dict__.iteritems():
          if getattr(pt.__dict__[key], 'atomic_number', -1) == type:
            type = key
            break
      self.add_atom(pos=line[1:4], type=type)

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
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}'.format(self)
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    self.hkl = array(value.split()[:3], dtype='int32')
  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}])'     \
           .format(self)

class Slabcut(Slabinfo):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl, isup, nl):
    super(Slabinfo, self).__init__(hkl)
    self.isup = isup
    """ Subsurface atomic label. """
    self.nl = nl
    """ Number of monolayers. """
  @property
  def raw(self):
    """ Input to SLABINFO """
    return '{0.hkl[0]} {0.hkl[1]} {0.hkl[2]}\n{0.isup} {0.nl}\n'.format(self)
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    line0, line1 = value.split('\n')[:2]
    self.hkl = array(line0.split()[:3], dtype='int32')
    self.isup, self.nl = [int(u) for u in line1.split()[:2]]

  def __repr__(self):
    return '{0.__class__.__name__}([{0.hkl[0]},  {0.hkl[1]}, {0.hkl[2]}], '    \
           '{0.isup}, {0.nl})'.format(self)

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
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    value = value.split('\n')
    n = int(value.pop(0).split()[0])
    for line in value[:n]: 
      line = line.split()
      type = int(line[0])
      pos = array(line[1:4], dtype='float64')
      self.add_atom(type=type, pos=pos)
    
   
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
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    from .. import periodic_table as pt
    value = value.split('\n')
    n = int(value.pop(0).split()[0])
    for line in value[:n]: 
      line = line.split()
      type = int(line[0])
      if type < 100: 
        for key, value in pt.__dict__.iteritems():
          if getattr(pt.__dict__[key], 'atomic_number', -1) == type:
            type = key
            break
      pos = array(line[1:4], dtype='float64')
      self.add_atom(type=type, pos=pos)

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
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    value = value.split('\n')
    n = int(value.pop(0).split())
    self.atoms = [u for v in value for u in v.split()]
    self.atoms = [int(u) for u in self.atoms[:n]]

class AffineTransform(Keyword):
  """ Affine transformation applied to the crystal or crystal fragment. """
  keyword = 'atomrot'
  """ CRYSTAL keyword. """
  def __init__( self, which=None, vectrans=None, origin=None, bondtrans=None,
                euler=None, bondrot=None, bondtoz=None ):
    """ Creates rotation. """
    self.which = which
    """ Selects what to apply transform to. 

        Should be either:

          - 'all' or None or an empty sequence, for the whole crystal.
          - a sequence of integers signifying the atoms to rotate
          - a single integer or a sequence of one integer, signifying that the
            molecule to which this atom belongs will be rotated.
    """
    from ..error import ValueError
    if [vectrans is None, origin is None, bondtrans is None].count(True) > 1:
      raise ValueError('More than one type of translation was selected.')
    if [euler is None, bondrot is None, bondtoz is None].count(True) > 1:
      raise ValueError('More than one type of rotation was selected.')

    self.vectrans  = vectrans
    self.bondtrans = bondtrans
    self.origin    = origin
    self.euler     = euler
    self.bondrot   = bondrot
    self.bondtoz   = bondtoz

  @property
  def vectrans(self):
    """ Selects translation by a vector. 

        Should be None (unselected) or a sequence of three real numbers
        defining the translation.
    """ 
    return self._vectrans
  @vectrans.setter
  def vectrans(self, value): 
    from numpy import array
    if value is None: self._vectrans = None; return
    self._vectrans = array(value, dtype='float64')
    self._origin, self._bondtrans = None, None
  @property
  def bondtrans(self):
    """ Selects translation by a vector. 

        Should be None (unselected) or a sequence of two labels(integers) and a
        real number. The first two indicate the axis and direction of the
        translation (first -> second), and the third its modulus in angstrom.

        .. _quantities: http://packages.python.org/quantities/index.html
    """ 
    return self._bondtrans
  @bondtrans.setter
  def bondtrans(self, value): 
    self._bondtrans = value
    if value is not None:
      self._origin, self._vectrans = None, None
  @property
  def origin(self):
    """ Shifts origin to specific atom. 
    
        Can be None or an atom label.
    """
    return self._origin
  @origin.setter
  def origin(self, value):
    if value is None: self._origin = None; return
    self._origin = int(value) 
    self._bondtrans, self._vectrans = None, None

  @property
  def euler(self):
    """ Defines a rotation using Euler matrices and a given atomic origin. 
    
        Should consist of three real numbers defining the Euler rotations (in
        degrees), and an atomic index defining the origin of the cartesian
        system. 
    """
    return self._euler
  @euler.setter
  def euler(self, value):
    self._euler = value
    if value is not None:
      self._bondrot, self._bondtoz = None, None
  @property
  def bondrot(self):
    """ Defines a rotation via two atoms forming the axis. 
 
        Should be two atomic indices followed by the rotation angle in degrees.
    """
    return self._bondrot
  @bondrot.setter
  def bondrot(self, value):
    self._bondrot = value
    if value is not None:
      self._euler, self._bondtoz = None, None
  @property
  def bondtoz(self):
    """ Defines a rotation via two atoms forming the axis. 
 
        Should be two atomic indices followed by the rotation angle in degrees.
    """
    return self._bondtoz
  @euler.setter
  def bondtoz(self, value):
    self._bondtoz = value
    if value is not None:
      self._euler, self._bondrot = None, None

  @property
  def raw(self):
    """ Creates raw input for crystal. """
    from quantities import angstrom, degree
    # first line
    if self.which is None: result = '0\n'
    elif not hasattr(self.which, '__iter__'): result = str(-self.which) + '\n'
    elif len(self.which) == 0: result = '0\n'
    elif len(self.which) == 1: result = str(-self.which[0]) + '\n'
    else:
      result = '{0}\n{1}\n'.format( len(self.which),                           \
                                    ' '.join(str(u) for u in self.which) )
    if self.origin is not None: result += str(self.origin) + ' '
    elif self.bondtrans is not None: result += '0 '
    elif self.vectrans is not None: result += '-1 '
    else: result += '999 '

    if self.bondtoz is not None: result += '1\n'
    elif self.bondrot is not None: result += '-1\n'
    elif self.euler is not None: result += '-2\n'
    else: result += '999\n'

    if self.bondtrans is not None:
      a, b, mod = int(self.bondtrans[0]), int(self.bondtrans[1]), self.bondtrans[2]
      if hasattr(mod, 'rescale'): mod = float(mod.rescale(angstrom).magnitude)
      result += '{0} {1} {2}\n'.format(a, b, mod)
    elif self.vectrans is not None:
      vec = self.vectrans
      if hasattr(vec, 'rescale'): vec = vec.rescale(angstrom).magnitude
      result += '{0[0]} {0[1]} {0[2]}\n'.format(vec)

    if self.bondtoz is not None: 
      a, b = int(self.bondtoz[0]), int(self.bondtoz[1])
      result += '{0} {1} 0\n'.format(a, b)
    elif self.bondrot is not None:
      a, b = int(self.bondrot[0]), int(self.bondrot[1])
      c = self.bondrot[2]
      if hasattr(c, 'rescale'): c = int(c.rescale(degree).magnitude + 0.01)
      result += '{0} {1} {2}\n'.format(a, b, c)
    elif self.euler is not None: 
      a, b, c = self.euler[:3]
      if hasattr(a, 'rescale'): a = int(a.rescale(degree).magnitude + 0.01)
      if hasattr(b, 'rescale'): b = int(b.rescale(degree).magnitude + 0.01)
      if hasattr(c, 'rescale'): c = int(c.rescale(degree).magnitude + 0.01)
      d = int(self.euler[4])
      result += '{0} {1} {2} {3}\n'.format(a, b, c, d)
    return result

  @raw.setter
  def raw(self, value):
    """ Read crystal input. """ 
    from numpy import array
    from quantities import angstrom, degree
     
    self.vectrans, self.bondtrans, self.origin = None, None, None
    self.euler, self.bondrot, self.bondtoz = None, None, None

    value = value.split()
    line = int(value.pop(0).split()[0])
    if line == 0: self.which = None
    elif line > 0: 
      self.which = []
      while len(self.which) < line:
        self.which.extend([u for u in value.pop(0).split()])
      self.which = [int(u) for u in self.which[:line]]
    else: self.which = line
    trans, rot = [int(u) for u in value.pop(0).split()[:2]]

    if trans > 0 and trans != 999: self.origin = trans
    elif trans == 0:
      line = value.pop(0).split()[:3]
      self.bondtrans = int(value[0]), int(value[1]), value[1] * angstrom
    else: 
      self.vectrans = array(value.pop(0).split()[:3], dtype='float64')         \
                      * angstrom        
    
    if rot < 0: 
      line = value.pop(0).split()[:4]
      self.euler = float(line[0]) * degree,                                    \
                   float(line[1]) * degree,                                    \
                   float(line[2]) * degree,                                    \
                   int(line[3])
    else:
      line = value.pop(0).split()[:4]
      a, b, alpha = int(line[0]), int(line[1]), int(line[2])
      if alpha == 0: self.bondtoz = a, b
      else: self.bondrot = a, b, alpha * degree

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
