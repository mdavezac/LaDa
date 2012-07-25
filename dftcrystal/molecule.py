__docformat__ = "restructuredtext en"
__all__ = ['Molecule']
from .input import ListBlock, Keyword, GeomKeyword
class Molecule(ListBlock):
  """ Molecule for the CRYSTAL code """
  keyword = 'molecule'
  def __init__(self, symmgroup=1, **kwargs):
    """ Creates a crystal following the CRYSTAL philosophy. """
    super(Molecule, self).__init__()

    self.symmgroup = symmgroup
    """ Symmetry group. """
    self.operators = kwargs.pop('operators', [])
    """ Functional operations on the structure. """
    self.atoms = []
    """ List of atoms. """
    for key, value in kwargs.iteritems(): setattr(self, key, value)

  def add_atom(self, *args, **kwargs):
    """ Adds an atom to the structure. """
    from ..crystal import Atom
    atom = Atom(*args, **kwargs)
    self.atoms.append(atom)
    return self

  def __repr_header__(self):
    """ Dumps representation to string. """
    # prints construction part.
    length = len('{0.__class__.__name__}('.format(self)) + 1
    args = [repr(self.symmgroup)]
    indent = ''.join([' ']*length) 
    for key, value in self.__dict__.iteritems():
      if key in ['operators', 'atoms', 'symmgroup'] : continue
      args.append('\\\n{2}{0}={1!r}'.format(key, value, indent))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __repr__(self):
    """ Dumps representation to string. """
    from inspect import getargspec
    result = self.__repr_header__()
    indent = ''.join([' '] * result.find('('))
    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent, dummy)
    # prints inner keywords/blocks
    for item in self: 
      hasindent = getargspec(item.__repr__).args
      hasindent = False if hasindent is None else 'indent' in hasindent
      result += '\\\n{0}.append({1!r})'.format(indent, item) if not hasindent  \
                else '\\\n{0}.append({1})'                                     \
                     .format(indent, item.__repr__(indent+'        ')) 
    return result

  @property 
  def raw(self):
    """ Raw input for structure. """
    from ..error import ValueError
    from ..periodic_table import find as find_specie

    
    # symmgroup 
    result = str(self.symmgroup).rstrip().lstrip() + '\n'

    # number of atoms + atoms
    result += '{0}\n'.format(len(self.atoms))
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
    if not hasattr(value, '__iter__'): value = value.split('\n')
    self.symmgroup = int(value.pop(0).split())

    n = int(value.pop(0).split()[0])
    self.atoms = []
    for line in value[:n]:
      line = line.split()
      type = int(line[0])
      if type < 100: type = find_specie(atomic_number=type).symbol
      self.add_atom(pos=[float(u) for u in line[1:4]], type=type)

  def read_input(self, tree, owner=None):
    """ Parses an input tree. """
    from parse import InputTree
    from . import registered
    self[:] = []
    self.raw = tree.raw
    do_breaksym = False
    has_breakkeep = False
    for key, value in tree:
      # parses sub-block.
      if isinstance(value, InputTree): continue
      # check for special keywords.
      if key.lower() == 'keepsymm':
        do_breaksym = False
        has_breakkeep = True
        continue
      if key.lower() == 'breaksym':
        do_breaksym = True
        has_breakkeep = True
        continue
      # finally, creates new object.
      if key.lower() in registered:
        newobject = registered[key.lower()]()
        if hasattr(newobject, 'breaksym'): newobject.breaksym = do_breaksym
      elif has_breakkeep: newobject = GeomKeyword(keyword=key, breaksym=do_breaksym)
      else: newobject = Keyword(keyword=key)
      if len(value) > 0: 
        try: newobject.raw = getattr(value, 'raw', value)
        except: pass
      self.append(newobject)
      do_breaksym = False
      has_breakkeep = False

  def append(self, keyword=None, raw=None, breaksym=None):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    from ..error import ValueError
    if isinstance(keyword, str):
      if breaksym is None: keyword = Keyword(keyword=keyword, raw=raw)
      else:
        breaksym = breaksym == True
        keyword = Keyword(keyword=keyword, raw=raw, breaksym=breaksym)
    elif not isinstance(keyword, Keyword):
      raise ValueError('Wrong argument to append.')
    list.append(self, keyword)
    return self

  def eval(self):
    """ Evaluates current structure. 
    
        Runs crystal to evaluate current structure.

        .. note:: Doesn't work yet.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. import crystal_program as program
    from .extract import Extract
    this = deepcopy(self)
    this.append('TESTGEOM')
    try:
      tmpdir = mkdtemp()

      with Changedir(tmpdir) as cwd:
        # writes input file
        with open('crystal.input', 'w') as file:
          file.write('{0}\n'.format(getattr(self, 'name', '')))
          file.write(this.print_input())
        # creates process and run it.
        if hasattr(program, '__call__'): program = program(self)
        process = Process( program, stdin='crystal.input',
                           outdir=tmpdir, stdout='crystal.out',
                           stderr='crystal.err', dompi=False )

        process.start()
        process.wait()

        return Extract().input_structure
    finally:
      try: rmtree(tmpdir)
      except: pass

  def copy(self):
    """ Returns deep copy of self. """
    from copy import deepcopy
    return deepcopy(self)
