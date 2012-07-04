__docformat__ = "restructuredtext en"
__all__ = ['BasisSet', 'Shell']
from quantities import UnitQuantity, angstrom
from .input import AttrBlock

crystal_bohr = UnitQuantity('crystal_bohr', 0.5291772083*angstrom, symbol='bohr')
""" Bhor radius as defined by CRYSTAL. """
crystal_invbohr2 = UnitQuantity( 'crystal_invbohr2',
                                 1e0/crystal_bohr/crystal_bohr, symbol='bohr^-2' )

class Shell(object):
  """ Defines a gaussian basis set for a specific orbital shell. """
  def __init__(self, type='s', charge=None, **kwargs):
    """ Creates a gaussian basis set. 
    
        :params str type:
          Denotes the orbital type of this shell.
        :param float charge:
          Denotes the formal electronic charge on this shell. 
          Defaults(None) to a full shell.
    """
    from operator import itemgetter
    super(Shell, self).__init__()
    self.type = type
    """ Orbital type of this shell """
    self.charge = charge
    """ Nominal charge of the shell """
    self.functions = []
    """ Radial functions description """
    # collect basis function keywords.
    keys = []
    for key, value in kwargs.iteritems():
      if key[0] != 'a':
        raise ValueError( 'Unexpected keyword {0!r} in Shell.__init__(...)'    \
                          .format(key) )
      try: int(key[1:])
      except:  
        raise ValueError( 'Unexpected keyword {0!r} in Shell.__init__(...)'    \
                          .format(key) )
      try: value = list(value)
      except: 
        raise ValueError( 'Unexpected keyword value for {0!r} '                \
                          'in Shell.__init__(...): {1!r}'.format(key, value) )
      if len(value) != (3 if self.type == 'sp' else 2):
        raise ValueError( 'Unexpected sequence length for {0!r} '              \
                          'in Shell.__init__(...): {1!r}'.format(key, value) )
      keys.append((key, int(key[1:]), value))
    for key, i, value in sorted(keys, key=itemgetter(1)):
      self.append(*value)
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
    if hasattr(alpha, 'rescale'): alpha = alpha.rescale(crystal_invbohr2)
    else: alpha = float(alpha) * crystal_invbohr2
    if self._type == 'sp': 
      if coef2 is None:
        raise TypeError('Expected second coeff for sp contraction.')
      self.functions.append([alpha, float(coef1), float(coef2)])
    else:
      if coef2 is not None:
        raise TypeError('Cannot use second coeff for in non-sp contraction.')
      self.functions.append([alpha, float(coef1)])
    return self

  @property
  def raw(self):
    """ Prints CRYSTAL input. """
    result = '0 {0.type_as_int} {1} {0.charge} 1\n'                            \
             .format(self, len(self.functions))
    for function in self.functions:
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbohr2).magnitude)
      if len(function) == 2:
        result += '{0:<20.12f}    {1: 20.12f}\n'.format(alpha, coef1)
      else:
        result += '{0:<20.12f}    {1:> 20.12f}  {2[2]:> 20.12f}\n'               \
                  .format(alpha, coef1, function)
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
    args = ['{0.type!r}'.format(self)]
    if self.charge != {'s': 2, 'p': 6, 'sp': 8, 'd': 10, 'f': 18}[self.type]:
      args.append('{0.charge!r}'.format(self))
    for i, function in enumerate(self.functions): 
      alpha, coef1 = function[0], function[1]
      if hasattr(alpha, 'rescale'):
        alpha = float(alpha.rescale(crystal_invbohr2).magnitude)
      vals = [alpha, coef1]
      if len(function) == 3: vals.append(function[2])
      args.append('a{0}={1!r}'.format(i, vals))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __len__(self): return len(self.functions)

class BasisSet(AttrBlock):
  """ Basis set block. """
  __ui_name__ = 'basis'
  """ Name used when printing user-friendly repr. """
  def __init__(self):
    """ Creates basis set block. """
    from .input import BoolKeyword, VariableListKeyword
    super(BasisSet, self).__init__()
    self._functions = {}
    """ Dictionary holding basis functions """
    self.atomsymm = BoolKeyword()
    """ Prints point symmetry for each atomic position """
    self.symmops  = BoolKeyword()
    """ Prints point symmetry operators """
    self.charged  = BoolKeyword()
    """ Allow charged system """
    self.noprint  = BoolKeyword()
    """ Disable basis set printing """
    self.paramprt = BoolKeyword()
    """ Print code dimensions parameters """
    self.ghosts   = VariableListKeyword(type=int)
    """ Remove atoms while keeping basis functions. """

  def _print_basis(self, structure=None):
    """ Prints only atoms in structure. """
    from ..error import KeyError
    from .. import periodic_table as pt
    # figure out set of species
    atoms = getattr(structure, 'atoms', structure)
    species = set([a.type for a in atoms]) if atoms is not None                \
              else self.iterkeys()
    # Now print them
    result = ''
    for key in species:
      if key not in self:
        raise KeyError('Could not find specie {0} in basis set.'.format(key))
      value = self[key]
      if len(value) == 0: continue
      if isinstance(key, str): key = getattr(pt, key).atomic_number
      result += '{0} {1}\n'.format(key, len(value))
      for function in value:
        result += function.raw.rstrip()
        if result[-1] != '\n': result += '\n'
      result = result.rstrip()
      if result[-1] != '\n': result += '\n'
    result += '99 0 end of basis functions per se\n'
    return result


  @property
  def raw(self):
    """ Raw basis set input. """
    return self._print_basis()
  @raw.setter
  def raw(self, value):
    """ Sets basis data from raw CRYSTAL input. """
    from ..error import NotImplementedError as NA, IOError
    # clean up all basis functions. 
    for key in self.iterkeys(): del self[key]
    lines = value.split('\n')
    while len(lines):
      line = lines.pop(0).split()
      type, nshells = int(line[0]), int(line[1])
      if type == 99: break
      type = self._specie(type)
      self._functions[type] = []
      for i in xrange(nshells):
        bt = int(lines[0].split()[0])
        if bt == 0:
          ncoefs = int(lines[0].split()[2])
          if len(lines) < ncoefs + 1: raise IOError('Unexpected end-of-file.')
          basis = Shell()
          basis.raw = '\n'.join(lines[:ncoefs+1])
          self._functions[type].append(basis)
          lines = lines[ncoefs+1:]
        else:
          raise NA('Basis type {0} hasnot yet been implemented.'.format(bt))


  def _specie(self, specie):
    """ Translate from specie to dictionary keyword. """
    from ..periodic_table import find as find_specie
    try: n = int(specie)
    except: pass
    else: specie = n
    if isinstance(specie, str): 
      specie = specie.rstrip().lstrip()
      try: specie = find_specie(name=specie)
      except: pass
      else: specie = specie.symbol
    elif specie < 100:
      try: specie = find_specie(atomic_number=specie)
      except: pass
      else: specie = specie.symbol
    return specie

  def __setitem__(self, specie, shell):
    """ Sets basis function for a particular specie. """
    specie = self._specie(specie)
    self._functions[specie] = shell
  def __getitem__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = self._specie(specie)
    return self._functions[specie]
  def __delitem__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = self._specie(specie)
    del self._functions[specie]
  def __contains__(self, specie):
    """ Gets basis-functions for a given specie. """
    specie = self._specie(specie)
    return specie in self._functions
  def __len__(self):
    """ Number of atomic species. """
    return len(self._functions)
  def __iter__(self):
    """ Iterates over species. """
    for key in self._functions: yield key
  def keys(self):
    """ List of species. """
    return self._functions.keys()
  def items(self):
    """ Tuples (specie, basis-functions) """
    return self._functions.items()
  def values(self):
    """ List of basis functions, """
    return self._functions.values()
  def iterkeys(self):
    """ Iterates over species. """
    return self._functions.iterkeys()
  def iteritems(self):
    """ Iterates over tuple (key, basis-functions). """
    return self._functions.iteritems()
  def itervalues(self):
    """ Iterates over basis-functions. """
    return self._functions.itervalues()

  def append(self, specie, shell):
    """ Adds shell to set. """
    specie = self._specie(specie)
    # now adds shell correctly.
    if specie in self._functions: self._functions[specie].append(shell)
    else: self._functions[specie] = [shell]
    return self

  def __repr__(self, defaults=False, name=None):
    """ Returns representation of this instance """
    from ..functools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def print_input(self, **kwargs):
    """ Dumps CRSTAL input to string. """
    post = super(BasisSet, self).print_input(**kwargs).rstrip()
    post = post[post.find('end of basis functions per se')+1:]
    post = post[post.find('\n')+1:] if post.find('\n') != -1 else ''
    result = self._print_basis(kwargs.get('structure', None)) + post
    if result[-1] != '\n': result += '\n'
    return result + 'END\n'

  def read_input(self, tree, owner=None):
    """ Parses an input tree. """
    self._functions = {}
    self._crysinput = self.__class__()._crysinput
    super(BasisSet, self).read_input(tree, owner=owner)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Prettier representation """
    from ..functools.uirepr import add_to_imports
    results = super(BasisSet, self).__ui_repr__(imports, name, defaults, ['raw'])
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    for key, value in self._functions.iteritems():
      newname = '{0}[{1!r}]'.format(name, key)
      results[newname] = '[{0}]'.format(', '.join(repr(u) for u in value))
      for u in value: add_to_imports(u, imports)
    return results
