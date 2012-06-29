__docformat__ = "restructuredtext en"
__all__ = ['Electronic']
from input import AttrBlock, TypedKeyword, BoolKeyword, Keyword

class Shrink(Keyword):
  """ k-point description -- SHRINK
  
      The IS (or IS1, IS2, IS3) and ISP keywords are mapped to
      :py:attr:`~Shrink.mp` and :py:att:`~Shrink.gallat`. 

      They can be used as follows::

        functional.shrink.mp = 5
        functional.shrink.gallat = 10

      This will map IS to 5 and ISP to 10.
      ISP automatically set to equal IS (or IS1) when :py:attr:`~Shrink.gallat`
      is set to None::

        functional.shrink.mp = 10
        functional.shrink.gallat = None

      This will print in the CRYSTAL input:

        | SHRINK
        | 5 5

      Finally, setting :py:attr:`~Shrink.mp` to a sequence of at most three
      integers will set IS to 0 and IS1, IS2 (defaults to 1), and IS3 (defaults
      to 1) to the relevant values::

        functional.shrink.mp = [5, 6]
        functional.shrink.gallat = None
        
      This will lead to the following input, where ISP defaulted automatically
      to IS1:

        | SHRINK
        | 0 5
        | 5 6 1
  """
  keyword = 'shrink'
  """ Crystal keyword. """
  def __init__(self, mp=1, gallat=None):
    """ Initializes Shrink keyword. 

        :param mp: 
           - if an integer, corresponds to IS variable of SHRINK.
           - if a sequence of at most three integers, corresponds to IS1, IS2,
             IS3. In that case, IS is of course set to 0. If there are only one
             or two numbers, then IS2 and/or IS3 is set to one. If None, the
             set to 1.
        :param int gallat:
           Corresponds to ISP variable of SHRINK. If None, then it is set equal
           to IS1
    """
    super(Shrink, self).__init__()
    self.mp = mp
    self.gallat = gallat
  @property
  def mp(self):
    """ Defines the Monkhorst-Pack mesh. """
    return self._mp
  @mp.setter
  def mp(self, value):
    if value is None: self._mp = 1
    elif hasattr(value, '__iter__'):
      mp = [max(int(u), 1) for u in value]
      if len(mp) > 3:
        raise ValueError('k-point mesh defined by at most three integers.')
      elif len(mp) == 0: mp = 1
      self._mp = mp
    else: self._mp = max(int(value), 1)
  @property
  def gallat(self): 
    """ Defines Gallat mesh for density matrix and Fermi energy. """
    return self._gallat
  @gallat.setter
  def gallat(self, value): 
    self._gallat = None if value is None else int(value)

  @property
  def raw(self):
    """ Prints raw CRYSTAL input. """
    if hasattr(self.mp, '__iter__'):
      ISP = self.mp[0] if self.gallat is None else self.gallat
      IS = ' '.join(str(u) for u in self.mp + [1]*(3-len(self.mp)))
      return '0 {0}\n{1}\n'.format(ISP, IS)
    else:
      ISP = self.mp if self.gallat is None else self.gallat
      return '{0} {1}\n'.format(self.mp, ISP)
  @raw.setter
  def raw(self, value):
    """ Sets from raw CRYSTAL input. """
    value = value.split()
    IS, ISP = int(value[0]), int(value[1])
    if IS == 0:
      self.mp = [int(u) for u in value[2:5]]
      if self.mp[-1] == 1: self.mp = self.mp[:-1]
      if self.mp[-1] == 1: self.mp = self.mp[:-1]
      if self.mp[-1] == 1: 
        self.mp = 1
        self.gallat = None if ISP == self.mp else ISP
      else: self.gallat = None if ISP == self.mp[0] else ISP
    else: 
      self.mp = IS
      self.gallat = None if ISP == self.mp else ISP
  def __repr__(self):
    """ Representation of this instance. """
    args = []
    if self.mp == 1:
      if self.gallat is not None:
        args.append('gallat={0.gallat!r}'.format(self))
    else:
      args.append(repr(self.mp))
      if self.gallat is not None: args.append('{0.gallat!r}'.format(self))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))


class Electronic(AttrBlock):
  """ DFT attribute block. """ 
  def __init__(self):
    """ Creates the scf attribute block. """
    from .hamiltonian import Dft
    super(Electronic, self).__init__()
    self.maxcycle = TypedKeyword(type=int)
    """ Maximum number of electronic minimization steps """
    self.tolinteg = TypedKeyword(type=[int]*5)
    """ Integration truncation criteria """
    self.toldep   = TypedKeyword(type=int)
    """ Density matrix convergence criteria """
    self.tolpseud = TypedKeyword(type=int)
    """ Pseudopotential truncation criteria """
    self.toldee   = TypedKeyword(type=int)
    """ Total energy convergence criteria """
    self.testpdim = BoolKeyword()
    """ Stop after processing input and performin symmetry analysis """
    self.test     = BoolKeyword()
    """ Stop after printing ressource requirement """
    self.symadapt = BoolKeyword()
    """ Symmetry adapted bloch wavefunctions """
    self.savewf   = BoolKeyword()
    """ Save wavefunctions to disk """
    self.shrink   = Shrink()
    """ k-point definition """
    self.fmixing  = TypedKeyword(type=int)
    """ Fock mixing during electronic minimiztion """
    self.dft      = Dft()
    """ Holds definition of the DFT functional itself """
