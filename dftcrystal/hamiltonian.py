__docformat__ = "restructuredtext en"
__all__ = ['Exchange', 'Correlation', 'B3LYP', 'B3PW', 'PBE0', 'SOGGAXC', 'Choice']
from . input import Keyword, AttrBlock, StrChoice

class Exchange(StrChoice):
  keyword = 'EXCHANGE'
  """ Corresponding CRYSTAL keyword. """
  def __init__(self, value='lda'):
    values =  ['becke', 'lda', 'pbe', 'pbesol', 'pwgga', 'sogga', 'vbh', 'wcgga']
    super(Exchange, self).__init__(values, value)
    self.excludes = ['b3lyp', 'b3pw', 'pbe0', 'soggaxc']
    """ keywords incompatible with this one. """
  def __set__(self, instance, value): 
    """ Makes sure that excluded keywords are set to False. """
    from ..error import ValueError
    if value is None:
      raise ValueError('Exchange and Correlations cannot be None.')
    super(Exchange, self).__set__(instance, value)
    for u in self.excludes: 
      if u in instance._crysinput: instance._crysinput[u] = False
class Correlation(StrChoice):
  keyword = 'CORRELAT'
  """ Corresponding CRYSTAL keyword. """
  def __init__(self, value='lda'):
    super(Correlation, self).__init__('pz')
    self.values = [ 'lyp', 'p86', 'pbe', 'pbesol', 'pwgga', 'pwlsd', 'pz',
                    'vbh', 'wl', 'vwn']

class Hybrid(Keyword):
  """ Adds hybrid exchange to DFT functional. """
  keyword = 'HYBRID'
  def __init__(self, value=0):
    super(Hybrid, self).__init__()
    self.raw = value
  @property
  def raw(self): return str(self._raw)
  @raw.setter
  def raw(self, value): self._raw = int(value)
  def __get__(self, instance, owner=None): return self._raw
  def __set__(self, instance, value):
    if value is None: self._raw = None
    self.raw = value
  def print_input(self, **kwargs):
    if self._raw is None: return None
    if self._raw == 0: return None
    return super(Hybrid, self).print_input(**kwargs)


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
  def raw(self): return str(self.exchange) + ' ' + str(self.correlation)
  @raw.setter
  def raw(self, value):
    value = value.split()
    self.exchange, self.correlation = float(value[0]), float(value[1])
  def disable(self): 
    """ Sets exchange and correlation to None. """
    self.exchange, self.correlation = None, None
  def print_input(self, **kwargs):
    if self.exchange is None or self.correlation is None: return None
    if abs(self.exchange) < 1e-8 and abs(self.correlation) < 1e-8: return None
    return super(Hybrid, self).print_input(**kwargs)


class Dft(AttrBlock):
  """ DFT attribute block. """ 
  keyword = 'DFT'
  def __init__(self):
    """ Creates the DFT attribute block. """
    super(Dft, self).__init__()
    self.exchange = Exchange() 
    """ Exchange functional. """
    self.correlation = Correlation() 
    """ Correlation functional. """
    self.hybrid = Hybrid()
    """ Amount of exchange to add to functional. """
    self.nonlocal = NonLocal()
    """ Non-local weights on exchange-correlation. """
    self.add_keyword('spin')
    self.spin = False
    """ If True, then perform spin-polarized calculation. """
  @property
  def b3lyp(self):
    """ Returns True if B3LYP functional is set. """
    if self.exchange != 'becke': return False
    if self.correlation != 'lyp': return False
    if self.hybrid is None: return False
    if int(self.hybrid) != 20: return False
    if self.nonlocal.exchange is None: return False
    if self.nonlocal.correlation is None: return False
    if abs(self.nonlocal.exchange - 0.9) > 1e-8: return False
    if abs(self.nonlocal.correlation - 0.81) > 1e-8: return False
    return True
  @b3lyp.setter
  def b3lyp(self, value):
    """ Sets functional to B3LYP if value is True. """
    if value == False: return
    self.exchange = 'becke'
    self.correlation = 'lyp'
    self.hybrid = 20
    self.nonlocal.exchange = 0.9
    self.nonlocal.correlation = 0.81
      
  @property
  def b3pw(self):
    """ Returns True if B3PW functional is set. """
    if self.exchange != 'becke': return False
    if self.correlation != 'pwgga': return False
    if self.hybrid is None: return False
    if int(self.hybrid) != 20: return False
    if self.nonlocal.exchange is None: return False
    if self.nonlocal.correlation is None: return False
    if abs(self.nonlocal.exchange - 0.9) > 1e-8: return False
    if abs(self.nonlocal.correlation - 0.81) > 1e-8: return False
    return True
  @b3pw.setter
  def b3pw(self, value):
    """ Sets functional to B3PW if value is True. """
    if value == False: return
    self.exchange = 'becke'
    self.correlation = 'pwgga'
    self.hybrid = 20
    self.nonlocal.exchange = 0.9
    self.nonlocal.correlation = 0.81
      
  @property
  def pbe0(self):
    """ Returns True if PBE0 functional is set. """
    if self.exchange != 'pbe': return False
    if self.correlation != 'pbe': return False
    if self.hybrid is not None: return False
    if self.hybrid != 0: return False
    if self.nonlocal.exchange is not None: return False
    if self.nonlocal.correlation is not None: return False
    if abs(self.nonlocal.exchange) > 1e-8: return False
    if abs(self.nonlocal.correlation) > 1e-8: return False
    return True
  @pbe0.setter
  def pbe0(self, value):
    """ Sets functional to PBE0 if value is True. """
    if value == False: return
    self.exchange = 'pbe'
    self.correlation = 'pbe'
    self.hybrid = 0
    self.nonlocal.exchange = 0
    self.nonlocal.correlation = 0
      
  @property
  def soggaxc(self):
    """ Returns True if SOGGAXC functional is set. """
    if self.exchange != 'sogga': return False
    if self.correlation != 'pbe': return False
    if self.hybrid is not None: return False
    if self.hybrid != 0: return False
    if self.nonlocal.exchange is not None: return False
    if self.nonlocal.correlation is not None: return False
    if abs(self.nonlocal.exchange) > 1e-8: return False
    if abs(self.nonlocal.correlation) > 1e-8: return False
    return True
  @soggaxc.setter
  def soggaxc(self, value):
    """ Sets functional to SOGGAXC if value is True. """
    if value == False: return
    self.exchange = 'sogga'
    self.correlation = 'pbe'
    self.hybrid._raw = 0
    self.nonlocal.exchange = 0
    self.nonlocal.correlation = 0

  def add_keyword(self, name, value=None):
    """ Adds input keyword. 

        Takes care of some global keywords and name changes.
    """
    # takes care of some global keywords
    if name in ['b3lyp', 'b3pw', 'pbe0', 'soggaxc']:
      setattr(self, name, value)
      return self
    # one case where name is transformed
    if name == 'correlat': name = 'correlation'
    # now can do normal add keyword.
    return super(Dft, self).add_keyword(name, value)
