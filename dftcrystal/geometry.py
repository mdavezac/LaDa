__docformat__ = "restructuredtext en"
__all__ = [ 'RemoveAtoms', 'Slabinfo', 'Slabcut', 'Slab',
            'DisplaceAtoms', 'InsertAtoms', 'ModifySymmetry',
            'AffineTransform', 'Marker' ]
from .input import Keyword, GeomKeyword

class RemoveAtoms(GeomKeyword):
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
    super(RemoveAtoms, self).__init__(**kwargs)
    self.labels = copy(args)
    """ Indices of atoms to remove. """

  @property
  def raw(self):
    """ Returns input string. """
    return str(len(self.labels)) + '\n'                                        \
           + ' '.join(str(u) for u in self.labels) + '\n'
  @raw.setter
  def raw(self, value):
    """ Read crystal input. """
    value = value.split('\n')
    n = int(value[0].split()[0])
    i = 1
    self.labels = []
    while len(self.labels) < n:
      self.labels.extend([u for u in value[i].split()])
      i += 1
    self.labels = [int(u) for u in self.labels[:n]]

  def __repr__(self):
    result = '{0.__class__.__name__}({1}, '                                    \
             .format( self, ', '.join(str(u) for u in self.labels))
    if self.breaksym == True: result += 'keepsym=False, '
    return result[:-2] + ')'
    
class Slabinfo(Keyword):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl=None):
    from numpy import array
    super(Slabinfo, self).__init__()
    self.hkl = hkl
    """ Surface normal. """
    if hkl is not None: self.hkl = array(hkl, dtype='int32')
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
  def print_input(self, **kwargs):
    """ Prints input. """
    if self.hkl is None: return None
    return super(Slabinfo, self).print_input(**kwargs)

class Slabcut(Slabinfo):
  """ Creates a slab from a 3d periodic system. """
  keyword = 'slabcut'
  """ CRYSTAL keyword. """
  def __init__(self, hkl=None, isup=None, nl=None, raw=None):
    super(Slabcut, self).__init__(hkl)
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
Slab = Slabcut
""" Alias for :py:class`~lada.dftcrystal.geometry.Slabcut`. """

class DisplaceAtoms(GeomKeyword):
  """ Displaces atoms. """
  keyword = 'atomdisp'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    """ Creates a displacement field. """
    super(DisplaceAtoms, self).__init__(**kwargs)
    self.atoms = []
    """ Atomic displacements. """

  def add_atom(self, *args, **kwargs):
    """ Adds a displacement to a given atom.
    
        At present, atom.type should be an index to an atom in the structure.
    """
    from ..crystal import Atom
    self.atoms.append(Atom(*args, **kwargs))
    return self

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
    result = '{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
      result += '{0.type} {0.pos[0]} {0.pos[1]} {0.pos[2]}\n'.format(atom)
    return result
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    from numpy import array
    self.atoms = []
    value = value.split('\n')
    n = int(value.pop(0).split()[0])
    for line in value[:n]: 
      line = line.split()
      type = int(line[0])
      pos = array(line[1:4], dtype='float64')
      self.add_atom(type=type, pos=pos)
    
   
class InsertAtoms(DisplaceAtoms):
  """ insert atom into structure. """
  keyword = 'atominse'
  """ CRYSTAL keyword. """
  def __init__(self, **kwargs):
    super(InsertAtoms, self).__init__(**kwargs)

  @property
  def raw(self):
    """ Raw input to CRYSTAL. """
    from .. import periodic_table as pt
    # number of atoms + atoms
    result = '{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
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
    self.atoms = []
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
    self.groups = []
    """ Atoms for which to modify symmetries. """
    for o in args:
      self.groups.append(list(o) if hasattr(o, '__iter__') else [o])

  @property
  def raw(self):
    """ Raw CRYSTAL input. """
    result = '{0}\n'.format(sum(len(o) for o in self.groups))
    for i, labels in enumerate(self.groups):
      result += ' '.join('{0} {1}'.format(u, i+1) for u in labels) + '\n'
    return result
  @raw.setter
  def raw(self, value):
    """ Reads input. """
    value = value.split()
    n = int(value.pop(0))
    d = {}
    for label, flag in zip(value[:2*n:2], value[1:2*n+1:2]):
      if flag in d: d[flag].append(int(label))
      else: d[flag] = [int(label)]
    self.groups = d.values()
  def __repr__(self):
    """ Representation of this instance. """
    return '{0.__class__.__name__}({1})'                                       \
           .format(self, ', '.join(repr(u) for u in self.groups))

class AffineTransform(Keyword):
  """ Affine transformation applied to the crystal or crystal fragment. """
  keyword = 'atomrot'
  """ CRYSTAL keyword. """
  def __init__( self, labels=None, vectrans=None, origin=None, bondtrans=None,
                euler=None, bondrot=None, bondtoz=None):
    """ Creates rotation. """
    self.labels = labels
    """ Selects what to apply transform to. 

        Should be either:

          - 'all' or None or an empty sequence, for the whole crystal.
          - a sequence of integers signifying the atoms to rotate
          - a single integer or a sequence of one integer, signifying that the
            molecule to which this atom belongs will be rotated.
    """
    from ..error import ValueError
    super(AffineTransform, self).__init__()

    if [vectrans is None, origin is None, bondtrans is None].count(False) > 1:
      raise ValueError('More than one type of translation was selected.')
    if [euler is None, bondrot is None, bondtoz is None].count(False) > 1:
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
    from quantities import angstrom
    if value is None: self._vectrans = None; return
    if hasattr(value, 'rescale'): self._vectrans = value 
    else: self._vectrans = array(value, dtype='float64') * angstrom
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
  @bondtoz.setter
  def bondtoz(self, value):
    self._bondtoz = value
    if value is not None:
      self._euler, self._bondrot = None, None

  @property
  def raw(self):
    """ Creates raw input for crystal. """
    from quantities import angstrom, degree
    # first line
    if self.labels is None: result = '0\n'
    elif not hasattr(self.labels, '__iter__'): result = str(-self.labels) + '\n'
    elif len(self.labels) == 0: result = '0\n'
    elif len(self.labels) == 1: result = str(-self.labels[0]) + '\n'
    else:
      result = '{0}\n{1}\n'.format( len(self.labels),                          \
                                    ' '.join(str(u) for u in self.labels) )
    if self.origin is not None: result += str(self.origin) + ' '
    elif self.bondtrans is not None: result += '0 '
    elif self.vectrans is not None: result += '-1 '
    else: result += '999 '

    if self.bondtoz is not None: result += '1\n'
    elif self.bondrot is not None: result += '1\n'
    elif self.euler is not None: result += '-1\n'
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
      d = int(self.euler[3])
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
    line = int(value.pop(0))
    if line == 0: self.labels = None
    elif line > 0: 
      self.labels = [int(u) for u in value[:line]]
      value = value[line:]
    else: self.labels = line
    trans, rot = [int(u) for u in value[:2]]
    value = value[2:]

    if trans > 0 and trans != 999: self.origin = trans
    elif trans == 0:
      self.bondtrans = int(value[0]), int(value[1]), float(value[2]) * angstrom
      value = value[3:]
    elif trans != 999: 
      self.vectrans = array(value[:3], dtype='float64') * angstrom        
      value = value[3:]
    
    if rot < 0: 
      self.euler = float(value[0]) * degree,                                   \
                   float(value[1]) * degree,                                   \
                   float(value[2]) * degree,                                   \
                   int(value[3])
    elif rot != 999:
      a, b, alpha = int(value[0]), int(value[1]), int(value[2])
      if alpha == 0: self.bondtoz = a, b
      else: self.bondrot = a, b, alpha * degree

  def __repr__(self): 
    """ Representation of AffineTransform. """
    args = []
    if self.labels is None: args.append('labels=None')
    elif isinstance(self.labels, int):
      args.append('labels={0.labels}'.format(self))
    elif len(self.labels) == 0: args.append('labels=None')
    elif len(self.labels) == 1:
      args.append('labels={0.labels[0]}'.format(self))
    else: args.append('labels={0.labels}'.format(self))
    if self.vectrans is not None:
      args.append('vectrans={0!r}'.format(self.vectrans))
    elif self.bondtrans is not None:
      args.append('bondtrans={0!r}'.format(self.bondtrans))
    elif self.origin is not None:
      args.append('origin={0!r}'.format(self.origin))
    if self.euler is not None:
      args.append('euler={0!r}'.format(self.euler))
    elif self.bondrot is not None:
      args.append('bondrot={0!r}'.format(self.bondrot))
    elif self.bondtoz is not None:
      args.append('bondtoz={0!r}'.format(self.bondtoz))
    return "{0.__class__.__name__}(".format(self) + ', '.join(args) + ')'

class Elastic(GeomKeyword):
  """ Elastic deformation of the lattice """
  keyword = 'elastic'
  """ CRYSTAL keyword """
  def __init__(self, matrix=None, is_epsilon=True, constvol=False, **kwargs):
    """ Creates cell-shape deformation. """
    super(Elastic, self).__init__(**kwargs)
    self.matrix = None
    self.is_epsilon = True
    self.const_volume = False
  @property
  def raw(self):
    if self.matrix is None: return ""
    type = 2 if self.is_epsilon else 1
    if not self.const_volume: type = -type
    result = str(type) + '\n'
    for i in xrange(3):
      result += ' '.join(str(self.matrix[i,j]) for j in xrange(3)) + '\n'
    return result
  @raw.setter
  def raw(self, value):
    from numpy import array
    value = value.split('\n')
    type = int(value[0].rstrip().lstrip())
    self.is_epsilon = abs(type) == 2
    self.const_volume = type > 0
    self.matrix = array([u.split()[:3] for u in value[1:4]], dtype='float64')

  def __repr__(self):
    """ Representation of this object. """
    args = []
    if self.matrix is not None:
      args.append(repr(self.matrix.tolist()))
    if not self.is_epsilon:
      args.append('is_epsilon=False' if len(args) > 0 else 'False')
    if self.const_volume:
      args.append('constvol=True' if len(args) > 0 else 'True')
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))


  def print_input(self, **kwargs):
    if self.matrix is None: return None
    return super(Elastic, self).print_input()


class Marker(Keyword):
  """ Places a marker in the execution list. """
  def __init__(self, name=None):
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
  def print_input(self, **kwargs): return None
