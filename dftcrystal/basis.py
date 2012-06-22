__docformat__ = "restructuredtext en"
__all__ = ['BasisSet', 'Shell']
from .input import Block
from quantities import UnitQuantity, angstrom

crystal_bhor = UnitQuantity('crystal_bhor', 0.5291772083*angstrom, symbol='bhor')
""" Bhor radius as defined by CRYSTAL. """
crystal_invbhor2 = UnitQuantity( 'crystal_invbhor2',
                                 1e0/crystal_bhor/crystal_bhor, symbol='bhor^-2' )

class Shell(object):
  """ Defines a gaussian basis set for a specific orbital shell. """
  def __init__(self, type='s', charge=None):
    """ Creates a gaussian basis set. 
    
        :params str type:
          Denotes the orbital type of this shell.
        :param float charge:
          Denotes the formal electronic charge on this shell. 
          Defaults(None) to a full shell.
    """
    super(Shell, self).__init__()
    self.type = type
    self.charge = charge
    self.functions = []

  @property
  def type(self):
    """ Orbital shell of this function. """
    return self._type
  @type.setter
  def type(self, value):
    from ..error import ValueError
    if isinstance(value, int):
      if value == 0:   value = 's'
      elif value == 1: value = 'sp'
      elif value == 2: value = 'p'
      elif value == 3: value = 'd'
      elif value == 4: value = 'f'
      else: value = int(value)
    elif value not in ['s', 'p', 'd', 'f', 'sp']:
      try: value = int(value)
      except: raise ValueError('Wrong input to type ({0}).'.format(value))
      else:
        self.type = value
        return 
    self._type = value
  @property
  def type_as_int(self):
    """ Orbital type as CRYSTAL input. """
    if self._type == 's': return 0
    elif self._type == 'sp': return 1
    elif self._type == 'p': return 2
    elif self._type == 'd': return 3
    elif self._type == 'f': return 4

  @property
  def charge(self):
    """ Formal electronic charge attributed to this shell. """
    return self._charge
  @charge.setter
  def charge(self, value):
    if value is None or value < 0:
      if self._type == 's': value = 2
      elif self._type == 'sp': value = 8
      elif self._type == 'p': value = 6
      elif self._type == 'd': value = 10
      elif self._type == 'f': value = 18
    self._charge = float(value)

  def append(self, alpha, coef1, coef2=None): 
    """ Appends a gaussian to the basis set. """
    from ..error import input, TypeError
    if len(self.functions) >= 10: 
      raise input('Cannot input further coefficients.') 
    if self._type == 'd' and len(self.functions) >= 6:
      raise input('Cannot input further coefficients.') 
    if hasattr(alpha, 'rescale'): alpha = alpha.rescale(crystal_invbhor2)
    else: alpha = float(alpha) * crystal_invbhor2
    if self._type == 'sp': 
      if coef2 is None:
        raise TypeError('Expected second coeff for sp contraction.')
      self.functions.append([alpha, coef1, coef2])
    else:
      if coef2 is not None:
        raise TypeError('Cannot use second coeff for in non-sp contraction.')
      self.functions.append([alpha, coef1])
    return self

  @property
  def raw(self):
    """ Prints CRYSTAL input. """
    result = '0 {0.type_as_int} {1} {0.charge} 1\n'                            \
             .format(self, len(self.functions))
    for function in self.functions:
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbhor2).magnitude)
      if len(function) == 2:
        result += '{0} {1}\n'.format(alpha, coef1)
      else: result += '{0} {1} {2[2]}\n'.format(alpha, coef1, function)
    return result

  @raw.setter
  def raw(self, value):
    """ Reads input lines. """
    lines = value.split('\n')
    dummy, self.type, n, self.charge, dummy = lines.pop(0).split()
    self.functions = []
    for i in xrange(int(n)):
      line = lines.pop(0).split()
      self.append(*line[:3 if self._type == 'sp' else 2]) 

  def __repr__(self, indent=''):
    """ Represention of this instance. """
    result = indent + '{0.__class__.__name__}({0.type!r}, {0.charge!r})'       \
                      .format(self)
    indent += '    '
    i = ' '.join(['']*result.find('('))
    for function in self.functions: 
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbhor2).magnitude)
      if len(function) == 2:
        args = '{0}, {1}'.format(alpha, coef1)
      else: args = '{0}, {1}, {2[2]}'.format(alpha, coef1, function)
      result += '\\\n{0}.append({1})'.format(i, args)
    return result

  def __len__(self): return len(self.functions)

class BasisSet(Block):
  """ Basis set block. """
  keyword = 'BASISSET'
  """ Fake CRYSTAL keyword. """
  def __init__(self):
    """ Creates basis set block. """
    super(BasisSet, self).__init__()
    self.basis = {}
    """ Mapping from species to atoms. """

  @property
  def raw(self):
    """ Raw basis set input. """
    from .. import periodic_table as pt
    result = ''
    for key, value in self.basis:
      if isinstance(key, str):
        type = getattr(pt, key, None)
        key = int(type) if type is None else type.atomic_number
      result += '{0} {1}\n{2}'.format(key, len(value), value.function.raw)
      result = result.rstrip()
      if result[-1] != '\n': result += '\n'
    result += '99 99\n'
    return result
  @raw.setter
  def raw(self, value):
    """ Sets basis data from raw CRYSTAL input. """
    from ..error import NotImplementedError as NA, IOError
    from ..periodic_table import find as find_specie
    self.basis = {}
    lines = value.split('\n')
    while len(lines):
      line = lines.pop(0).split()
      type, nshells = int(line[0]), int(line[1])
      if type == 99 or nshells == 99: break
      if type < 100: 
        try: type = find_specie(atomic_symbol=type).symbol
        except: pass
      self.basis[type] = []
      for i in xrange(nshells):
        bt = int(lines[0].split()[0])
        if bt == 0:
          ncoefs = int(lines[0].split()[2])
          if len(lines) < ncoefs + 1: raise IOError('Unexpected end-of-file.')
          basis = Shell()
          basis.raw = '\n'.join(lines[:ncoefs+1])
          self.basis[type].append(basis)
          lines = lines[ncoefs+1:]
        else:
          raise NA('Basis type {0} hasnot yet been implemented.'.format(bt))
  def add_shell(self, name, shell):
    """ Adds shell to set. """
    from ..periodict_table import find as find_specie

    # figure out what name is.
    try: n = int(name)
    except: pass
    else: name = n
    if isinstance(name, str): 
      name = name.rstrip().lstrip()
      try: specie = find_specie(name=name)
      except: pass
      else: name = specie.symbol
    elif name < 100:
      try: specie = find_specie(atomic_number=name)
      except: pass
      else: name = specie.symbol

    # now adds shell correctly.
    if name not in self.basis: self.basis[name] = [shell]
    else: self.basis[name].append(shell)
    
  def __repr__(self, indent=''):
    """ Representation of this instance. """
    result = indent + '{0.__class__.__name__}()'.format(self)
    indent += '    '
    for key, value in self.basis.iteritems():
      for func in value: 
        result += '\\\n{0}.add_shell({1!r}, {2!r}})'.format(indent, key, func)
    return result

  def print_input(self, **kwargs):
    """ Dumps CRSTAL input to string. """
    result = super(BasisSet, self).print_input(**kwargs)
    return result[result.find('\n')+1:]

  def read_input(self, tree):
    """ Parses an input tree. """
    from input import InputTree, Keyword
    from ..error import ValueError
    from . import registered
    self[:] = []
    self.raw = tree.raw
    for key, value in tree:
      # parses sub-block.
      if isinstance(value, InputTree):
        raise ValueError( '{0} is not allowed in the basis set definition.'    \
                          .format(key) )
      # creates new object.
      newobject = registered.get(key.lower(), Keyword(keyword=key))
      if len(value) > 0: 
        try: newobject.raw = getattr(value, 'raw', value)
        except: pass
      self.append(newobject)
