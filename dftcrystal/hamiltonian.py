__docformat__ = "restructuredtext en"
__all__ = ['Exchange', 'Correlation', 'NonLocal', 'Dft', 'GlobalExc', 'Hybrid']
from . input import Keyword, AttrBlock, Choice, TypedKeyword, BoolKeyword

class Exchange(Choice):
  keyword = 'exchange'
  """ Corresponding CRYSTAL keyword. """
  values = ['becke', 'lda', 'pbe', 'pbesol', 'pwgga', 'sogga', 'vbh', 'wcgga']
  """" Set of values from which to choose. """
  def __init__(self, value=None): 
    super(Exchange, self).__init__(value)
    self.excludes = ['b3lyp', 'b3pw', 'pbe0', 'soggaxc']
    """ keywords incompatible with this one. """
  def __set__(self, instance, value): 
    """ Makes sure that excluded keywords are set to False. """
    from ..error import ValueError
    if value is None:
      raise ValueError('Exchange and Correlations cannot be None.')
    super(Exchange, self).__set__(instance, value)
    for u in self.excludes: 
      if getattr(instance, u, False): setattr(instance, u, False)
  def print_input(self, **kwargs):
    """ Print EXCHANGE keyword. 

        Also checks whether one of global keywords is set, in which case does
        not print.
    """
    if self._value is None: return None
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    return '{0}\n{1}\n'.format(self.keyword.upper(), self.value.upper())

class Correlation(Exchange):
  keyword = 'correlat'
  """ Corresponding CRYSTAL keyword. """
  values = [ 'lyp', 'p86', 'pbe', 'pbesol', 'pwgga', 'pwlsd', 'pz',
             'vbh', 'wl', 'vwn']
  """ Set of values from which to choose. """
  def __init__(self, value=None):
    super(Correlation, self).__init__(value=value)

class GlobalExc(Keyword):
  """ Takes care of global keywords for exchange-correlation. """
  def __init__(self, keyword, values):
    super(GlobalExc, self).__init__(keyword=keyword)
    self.values = values
    """ Values against which to test. """

  def __get__(self, instance, owner=None):
    """ True if global keyword is set. """
    def test(a, b):
      if a is None:
        if b is not None: return False
      else: 
        if b is None: return False
        if abs(a-b) > 1e-8: return False
      return True
    return instance.exchange == self.values[0]                                 \
           and instance.correlat == self.values[1]                             \
           and test(self.values[2], instance.hybrid)                           \
           and test(self.values[3], instance.nonlocal.exchange)                \
           and test(self.values[4], instance.nonlocal.correlation)
  def __set__(self, inst, value):
    """ Sets to global if value is True. 
    
        If False, ignores.
    """
    if value == True:
      inst.exchange, inst.correlat, inst.hybrid,                               \
        inst.nonlocal.exchange, inst.nonlocal.correlation = self.values
  def print_input(self, **kwargs):
    if self.__get__(kwargs['crystal'].dft): 
      return self.keyword.upper() + '\n'
  def read_input(self, tree, owner=None):
    """ True if ever read. """
    self.__set__(owner, True)

class NonLocal(Keyword):
  """ Non-local weight parameters. """
  keyword = 'NONLOCAL'
  def __init__(self, xc=None, corr=None):
    super(NonLocal, self).__init__()
    self.exchange = None
    """ Exchange weight. """
    self.correlation = None
    """ Correlation weight. """
  @property
  def raw(self):
    x = 0 if self.exchange is None else self.exchange
    c = 0 if self.correlation is None else self.correlation
    return '{0} {1}\n'.format(x, c)
  @raw.setter
  def raw(self, value):
    value = value.split()
    self.exchange, self.correlation = float(value[0]), float(value[1])
  def disable(self): 
    """ Sets exchange and correlation to None. """
    self.exchange, self.correlation = None, None
  def print_input(self, **kwargs):
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    x = 0 if self.exchange is None else self.exchange
    c = 0 if self.correlation is None else self.correlation
    if abs(x) < 1e-8 and abs(c) < 1e-8: return None
    return super(NonLocal, self).print_input(**kwargs)
  def __getitem__(self, index):
    """ list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: return self.exchange
    elif index == 1 or index == -1: return self.correlation
    raise IndexError('nonlocal can be indexed with 0, 1, or -1 only.')
  def __setitem__(self, index, value):
    """ sets as list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: self.exchange = value
    elif index == 1 or index == -1: self.correlation = value
    else: raise IndexError('nonlocal can be indexed with 0, 1, or -1 only.')
  def __len__(self): return 2

class Hybrid(TypedKeyword):
  keyword = 'hybrid'
  """ Crystal input keyword. """
  type = float
  """ Keyword must conform to this type, if not None. """
  def __init__(self, value=None):
    super(Hybrid, self).__init__(value=value)
  def print_input(self, **kwargs):
    inst = kwargs['crystal'].dft
    if inst.pbe0 or inst.b3lyp or inst.b3pw or inst.soggaxc: return None
    return super(Hybrid, self).print_input(**kwargs)



class Dft(AttrBlock):
  """ DFT attribute block. """ 
  keyword = 'DFT'
  def __init__(self):
    """ Creates the DFT attribute block. """
    super(Dft, self).__init__()
    self.exchange = Exchange() 
    """ Exchange functional. """
    self.correlat = Correlation() 
    """ Correlation functional. """
    self.hybrid   = Hybrid()
    """ Amount of exchange to add to functional. """
    self.nonlocal = NonLocal()
    """ Non-local weights on exchange-correlation. """
    self.spin     = BoolKeyword(keyword='spin')
    """ If True, then perform spin-polarized calculation. """
    self.b3lyp    = GlobalExc('b3lyp', ['becke', 'lyp', 20, 0.9, 0.81])
    """ B3LYP global keyword. """
    self.b3pw     = GlobalExc('b3pw', ['becke', 'pwgga', 20, 0.9, 0.81])
    """ B3PW global keyword. """
    self.pbe0     = GlobalExc('pbe0', ['pbe', 'pbe', None, None, None])
    """ B3PW global keyword. """
    self.soggaxc  = GlobalExc('soggaxc', ['sogga', 'pbe', None, None, None])
    """ B3PW global keyword. """

  def __ui_repr__( self, imports, name=None,
                   rmdefaults=False, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    exclude = [] if exclude is None else list(exclude)
    if not self.b3lyp: exclude.append('b3lyp')
    if not self.b3pw: exclude.append('b3pw')
    if not self.pbe0: exclude.append('pbe0')
    if not self.soggaxc: exclude.append('soggaxc')
    if self.b3lyp or self.b3pw or self.pbe0 or self.soggaxc: 
      exclude.extend(['exchange', 'correlat', 'hybrid', 'nonlocal'])
    return super(Dft, self).__ui_repr__(imports, name, rmdefaults, defaults, exclude)

