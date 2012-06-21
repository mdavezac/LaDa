__docformat__ = "restructuredtext en"
__all__ = ['BasisSet', 'Shell']
from .input import Block

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
    self.type = type
    self.charge = charge
    self.functions = []

  @property
  def type(self):
    """ Orbital shell of this function. """
    return self._type
  @type.setter
  def type(self, value):
    if isinstance(type, int):
      if type == 0: type = 's'
      elif type == 1: type = 'sp'
      elif type == 2: type = 'p'
      elif type == 3: type = 'd'
      elif type == 4: type = 'f'
      else: type = int(type)
    self._type = type
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
    self.value = float(value)

  def append(self, alpha, coef1, coef2=None): 
    """ Appends a gaussian to the basis set. """
    from ..error import input
    from ..quantities import angstrom
    bhor = 0.5291772083*angstrom # crystal conversion factor
    if len(self.functions) >= 10: 
      raise input('Cannot input further coefficients.') 
    if self._type == 'd' and len(self.functions) >= 6:
      raise input('Cannot input further coefficients.') 
    if hasattr(alpha, 'rescale'): alpha = alpha.rescale(1./bhor/bhor)
    else: alpha *= 1e0/bhor/bhor
    if self._type == 'sp': 
      self.functions.append([alpha, coef1, coef2])
    else: self.functions.append([alpha, coef1])

  @property
  def raw(self):
    """ Prints CRYSTAL input. """
    from ..quantities import angstrom
    bhor = 0.5291772083*angstrom # crystal conversion factor
    invbhor = 1e0/bhor/bhor
    result = '0 {0.type_as_int} {1} {0.charge} 1\n'                            \
             .format(self, len(self.functions))
    for function in self.functions:
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(invbhor).magnitude)
      if len(function) == 2:
        result += '{0} {1}\n'.format(alpha, coef1)
      else: result += '{0} {1} {2[2]}\n'.format(alpha, coef1, function)

  def read_input(self, lines): 
    """ Reads and consumes input lines. """
    dummy, self.type, n, self.charge, dummy = lines.pop(0).split()
    self.functions = []
    for i in xrange(n):
      line = lines.pop(0).split()
      self.append(*line[:3 if self._type == 'sp' else 2]) 

  def __repr__(self, indent=''):
    """ Represention of this instance. """
    from ..quantities import angstrom
    bhor = 0.5291772083*angstrom # crystal conversion factor
    invbhor = 1e0/bhor/bhor
    result = indent + '{0.__class__.__name__}({0.type!r}, {0.charge!r})'       \
                      .format(self)
    indent += '    '
    i = result.find('(')
    for function in self.functions: 
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(invbhor).magnitude)
      if len(function) == 2:
        args = '{0}, {1}\n'.format(alpha, coef1)
      else: args = '{0}, {1}, {2[2]}\n'.format(alpha, coef1, function)
      result += '\\\n{0}.append({1})'.format(i, args)

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
    """ Prints raw basis set input. """
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
