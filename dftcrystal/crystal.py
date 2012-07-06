__docformat__ = "restructuredtext en"
__all__ = ['Crystal']
from .molecule import Molecule
class Crystal(Molecule):
  """ CRYSTAL-wise structure, e.g. functional approach. """
  keyword = 'crystal'
  """ CRYSTAL keyword. """
  def __init__(self, symmgroup=1, *args, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    shift = kwargs.pop('shift', 0)
    super(Crystal, self).__init__(symmgroup, **kwargs)

    self.params = list(args)
    """ Lattice parameters. """
    self.shift = shift

  @property
  def shift(self):
    """ Shift of the origin of the cell,
    
         - if 0 or second, second setting given international tables
           (equivalent to IFSO == 0).
         - if 1 or first, second setting given international tables
           (equivalent to IFSO == 1).
         - if a sequence of three, then an actual shift in reciprocal space * 24.
    """
    return self._shift
  @shift.setter
  def shift(self, value):
    from numpy import array
    from ..error import ValueError
    if hasattr(value, '__iter__'): 
      self._shift = array(value, dtype='float64')
    elif value == 0 or value == 'second': self._shift = 0
    elif value == 1 or value == 'first': self._shift = 1
    else: raise ValueError('Unknown input value for shift.')

  def __repr_header__(self):
    """ Dumps representation to string. """
    # prints construction part.
    length = len('{0.__class__.__name__}('.format(self)) + 1
    args = [repr(self.symmgroup)]
    for o in self.params: args.append('{0!r}'.format(o))
    indent = ' '.join('' for i in xrange(length))
    for key, value in self.__dict__.iteritems():
      if key in ['operators', 'atoms', 'symmgroup', 'params'] : continue
      if key == '_shift': key = 'shift'
      args.append('\\\n{2}{0}={1!r}'.format(key, value, indent))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  @property 
  def raw(self):
    """ Raw input for structure. """
    from ..physics.spacegroup import space_group as sg
    from ..error import ValueError
    from ..periodic_table import find as find_specie

    if hasattr(self.shift, '__iter__'): hasshift = 2
    else:                               hasshift = self.shift

    try: symmgroup = int(self.symmgroup)
    except: 
      if self.symmgroup not in sg:
        raise ValueError('Unknown space group {0.symmgroup}'.format(self))
      symmgroup = sg[self.symmgroup].index
    
    # control flags + symmgroup 
    result = '0 {0} {1}\n{2}\n'.format( getattr(self, 'ifhr', 0),              \
                                        hasshift, symmgroup )
    # shift if needed, watch out for weird ass input.
    if hasshift == 2:
      shift = self.shift * 24.0 + 0.1
      shift = [int(u) for u in shift]
      result += '{0[0]} {0[1]} {0[2]}\n'.format(shift)

    # lattice parameters.
    for o in self.params: result += '{0} '.format(o)

    # number of atoms + atoms
    result += '\n{0}\n'.format(len(self.atoms))
    for atom in self.atoms:
      try: n = find_specie(name=atom.type).atomic_number
      except:
        try: n = int(atom.type)
        except: 
          raise ValueError( 'Could not make sense of atomic type {0.type}.'    \
                            .format(atom) )
      result += '{0} {1.pos[0]} {1.pos[1]} {1.pos[2]}\n'.format(n, atom)
    return result

  @raw.setter
  def raw(self, value):
    """ Reads crystal input. """
    from ..periodic_table import find as find_specie
    from numpy import array
    if not hasattr(value, '__iter__'): value = value.split('\n')
    data = value.pop(0).split()
    self.symmgroup, self.ifhr, shift = (int(u) for u in data[:3]) 
    data = value.pop(0)
    if self.symmgroup == 0: self.symmgroup = int(data.split()[0])
    else: self.symmgroup = data

    if shift > 1: 
      self.shift = array(value.pop(0).split(), dtype='float64') / 24.0
    else: self.shift = shift

    self.params = [float(u) for u in value.pop(0).split()]
    n = int(value.pop(0).split()[0])
    self.atoms = []
    for line in value[:n]:
      line = line.split()
      type = int(line[0])
      if type < 100: type = find_specie(atomic_number=type).symbol
      self.add_atom(pos=[float(u) for u in line[1:4]], type=type)
