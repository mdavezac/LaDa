__docformat__ = "restructuredtext en"
__all__ = ['External']
from .molecule import Molecule
class External(Molecule):
  """ Functional CRYSTAL-wise structure, starting from LaDa structure.

      Provides a mixture of CRYSTAL_ and LaDa-style structures. The starting
      structure is a :py:class:`~lada.crystal.cppwrappers.Structure` instance
      which is inputed into CRYSTAL_ via the EXTERNAL keyword. This instance
      also accepts functional modifications which will act upon the initial
      structure, much as :py:class:`~lada.dftcrystal.crystal.Crystal`
      instances.
  """
  keyword = 'external'
  """ CRYSTAL keyword. """
  def __init__(self, *args, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    from ..crystal import Structure
    self.__dict__['initial'] = None
    super(External, self).__init__()
    del self.symmgroup
    del self.__dict__['atoms']

    self.initial = kwargs.pop('copy', None)
    """ Initial structure. """
    if self.initial is None: self.initial = Structure(*args, **kwargs)
    else: self.initial = self.initial.copy()

  def add_atoms(self, *args, **kwargs): 
    """ Adds atom to initial structure. """
    self.initial.add_atom(*args, **kwargs)
    return self
  @property
  def atoms(self):
    """ Alias to the initial structure. """
    return self.initial
  @atoms.setter
  def atoms(self, value):
    """ Alias to the initial structure. """
    from ..error import AttributeError
    if self.initial is not None: 
      raise AttributeError('Cannot set atoms in External instance.')
    self.__dict__['atoms'] = value

  @property
  def cell(self):
    """ Cell-vectors of the initial instance. """
    return self.initial.cell
  @cell.setter
  def cell(self, value):
    self.initial.cell = value

  def __repr_header__(self):
    """ Dumps representation to string. """
    # prints construction part.
    args = ""
    if self.initial is not None:
      copy = self.initial.copy()
      copy.clear()
      args = repr(copy)
      args = args[args.find('(')+1:args.rfind(')')]
    return '{0.__class__.__name__}({1})'.format(self, args)

  @property 
  def raw(self):
    """ Raw input for structure. """
    return ""
  @raw.setter 
  def raw(self, value): pass

  def output_map(self, **kwargs):
    """ Writes external file and calls base class. """
    from ..crystal import write
    if kwargs.get('filework', False) == True:  write.crystal(self.initial)
    return super(External, self).output_map(**kwargs)
