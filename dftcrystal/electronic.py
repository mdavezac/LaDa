__docformat__ = "restructuredtext en"
__all__ = ['Electronic']
from .input import AttrBlock, BoolKeyword
from ..tools.input import TypedKeyword, BaseKeyword
from quantities import UnitQuantity, hartree

class Shrink(BaseKeyword):
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


      Another option it to set :py:attr:`~Shrink.mp` and
      :py:attr:`~Shrink.gallat` directly::

        functional.shrink = 8, 8

      would result in the following ouput

        | SHRINK
        | 8 8

      :py:attr:`~Shrink.mp` will be set to the first item. and
      :py:attr:`~Shrink.gallat` to the second. There should always be two
      items, unless setting to None.
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

  def __set__(self, instance, value):
    """ Sets the keyword more easily. """
    if value is None: self.mp = None; return
    self.mp = value[0]
    self.gallat = value[1]

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

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    from ..tools.uirepr import add_to_imports
    if defaults is not None: 
      if type(defaults) is not type(self): 
        add_to_imports(self)
        return {name: repr(self)}
      if self.gallat == defaults.gallat and self.mp == defaults.mp: 
        return {}
      return {name: '{0.mp!r}, {0.gallat!r}'.format(self)}
    elif name is None:
      add_to_imports(self, imports)
      return {None: 'shrink = {0!r}'.format(self)}
    add_to_imports(self, imports)
    return {name: self.__repr__()}

  def output_map(self, **kwargs):
    """ Prints SHRINK keyword """
    from .molecule import Molecule
    # The 'get' piece is to make testing a bit easier
    if type(kwargs.get('structure', None)) is Molecule: return None
    return super(Shrink, self).output_map(**kwargs)


class LevShift(BaseKeyword):
  """ Implements LevShift keyword. """
  keyword = 'levshift'
  """ CRYSTAL input keyword. """
  lock = BoolKeyword(value=True)
  """ If True, maintains shift after diagonalization.
  
      Defaults to True.
  """
  units = UnitQuantity('decihartree', 0.1*hartree)
  """ LEVSHIFT takes 0.1 * hartrees. """
  def __init__(self, shift=None, lock=None):
    """ Creates the LEVSHIFT keyword. """
    super(LevShift, self).__init__()
    self.shift = shift
    self.lock = lock
  @property
  def shift(self): 
    """ Value to which this keyword is set. """
    return self._shift
  @shift.setter
  def shift(self, shift):
    """ Artificial shift between conduction and valence bands. 
  
        Can be a regular floating point, in which case the units are 0.1 *
        hartree, or a scalar signed by an unit with dimensionality of an energy.
        In the latter case, the value is rescaled to 0.1 H.
    """
    from ..error import ValueError
    if shift is None: self._shift = None; return
    if hasattr(shift, 'rescale'): shift = shift.rescale(self.units)
    else: shift = shift * self.units
    if len(shift.shape) != 0:
      raise ValueError('shift should be a scalar.')
    self._shift = shift
  @property
  def lock(self): 
    """ Whether shift is kept after diagaonalization. """
    return self._lock
  @lock.setter
  def lock(self, lock):
    self._lock = None if lock is None else (lock == True)
  def __set__(self, instance, value):
    """ Sets the value of this instance. """
    from ..error import ValueError
    if not hasattr(value, '__getitem__'): 
      raise ValueError('Incorrect input to levshift: {0}.'.format(value))
    self.shift = value[0]
    self.lock = value[1]
  def __getitem__(self, index):
    """ list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: return self.shift
    elif index == 1 or index == -1: return self.lock
    raise IndexError('Levshift can be indexed with 0, 1, or -1 only.')
  def __setitem__(self, index, value):
    """ sets as list [self.shift, self.lock] """
    from ..error import IndexError
    if index == 0: self.shift = value
    elif index == 1 or index == -1: self.lock = value
    else: raise IndexError('Levshift can be indexed with 0, 1, or -1 only.')
  def __len__(self): return 2
  @property
  def raw(self):
    """ CRYSTAL input as a string """
    if self.shift is None or self.lock is None: return ''
    return '{0} {1}'.format(int(float(self.shift)+0.01), 1 if self.lock else 0)
  @raw.setter
  def raw(self, value):
    """ Sets instance from raw data """
    self.shift = int(value.split()[0])
    self.lock = int(value.split()[1]) != 0
  def output_map(self, **kwargs):
    if self.shift is None or self.lock is None: return None
    return super(LevShift, self).output_map(**kwargs)
  def __repr__(self):
    args = []
    if self.shift is not None: args.append(str(float(self.shift)))
    if self.lock is not None:
      args.append( 'lock={0}'.format(self.lock) if len(args) == 0              \
                   else str(self.lock) )
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

class GuessP(BoolKeyword):
  """ Reads density matrix from disk.

      If True *and* restart_ is not None, then copies crystal.f9 to
      fort.20 in the working directory and adds GUESSP keyword to the input.

      If True but restart_ is None or the file crystal.f9 does not exist, then
      does nothing. This is not an error, however.

      If False or None, does nothing. Since GuessP can lead to bizarre
      problems, it is *False* by default.

      .. _restart: :py:attr:`~lada.dftcrystal.functional.restart` 
  """ 
  keyword = 'guessp'
  def __init__(self, value=True):
    super(GuessP, self).__init__(value=value)
  def output_map(self, **kwargs):
    from os.path import exists, join
    from ..misc import copyfile
    if self.value is None or self.value == False: return None
    if kwargs['crystal'].restart is None: return None
    path = join(kwargs['crystal'].restart.directory, 'crystal.f9')
    if not exists(path): return None
    copyfile(path, join(kwargs['workdir'], 'fort.20'), nothrow='same')
    return super(GuessP, self).output_map(**kwargs)


class Electronic(AttrBlock):
  """ DFT attribute block. """ 
  __ui_name__ = 'electronics'
  """ Name used when printing this instance with ui_repr """
  def __init__(self):
    """ Creates the scf attribute block. """
    from ..tools.input import ChoiceKeyword
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
    self.levshift = LevShift()
    """ Fock mixing during electronic minimiztion """
    self.ppan     = BoolKeyword()
    """ Mulliken population analysis """
    self.biposize = TypedKeyword(type=int)
    """ Size of buffer for Coulomb integrals bipolar expansions """
    self.exchsize = TypedKeyword(type=int)
    """ Size of buffer for exchange integrals bipolar expansions """
    self.scfdir   = BoolKeyword()
    """ Whether to reevaluate integrals at each electronic step """
    self.poleordr = ChoiceKeyword(values=range(0, 7))
    """ Coulomb intergrals pole truncation """
    self.guessp   = GuessP(value=False)
    """ Reads density matrix from disk.
    
        If True *and* restart_ is not None, then copies crystal.f9 to
        fort.20 in the working directory and adds GUESSP keyword to the input.
    
        If True but restart_ is None or the file crystal.f9 does not exist, then
        does nothing. This is not an error, however.
    
	If False or None, does nothing. Since GuessP can lead to bizarre
        problems, it is *False* by default.
    
        .. _restart: :py:attr:`~lada.dftcrystal.functional.restart` 
    """ 
    self.dft      = Dft()
    """ Holds definition of the DFT functional itself """
    self.nofmwf   = BoolKeyword()
    """ Whether to print formatted wavefunctions to disk. """
    self.nobipola = BoolKeyword()
    """ Whether to compute bielectronic integrals exactly. """
    self.nobipcou = BoolKeyword()
    """ Whether to compute bielectronic Coulomb integrals exactly. """
    self.nobipexc = BoolKeyword()
    """ Whether to compute bielectronic exchange integrals exactly. """
    self.nomondir = BoolKeyword()
    """ Whether store monoelectronic integrals to disk. 
    
        If True, the monoelctronic integrals are computed once at the start of
        the SCF calculation and re-read from disk at each geometric
        optimization step.
    """
    self.nosymada = BoolKeyword()
    """ Whether to not use symmetry adapted functions. """
    self.savewf   = BoolKeyword()
    """ Whether to save wavefunctions at each step. """
    self.mpp      = BoolKeyword(value=False)
    """ Whether to use MPP or Pcrystal when running mpi. """
