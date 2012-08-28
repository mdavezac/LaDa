__docformat__ = "restructuredtext en"
__all__ = ['Molecule']
from .input import GeomKeyword
from ..tools.input import ListBlock, BaseKeyword
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
    """ Dumps representation header to string. 
    
        Mostly useful to limit how much should be rewritten of repr in derived
        classes.
    """
    # prints construction part.
    length = len('{0.__class__.__name__}('.format(self)) + 1
    args = [repr(self.symmgroup)]
    indent = ''.join([' ']*length) 
    for key, value in self.__dict__.iteritems():
      if key in ['operators', 'atoms', 'symmgroup'] : continue
      args.append('\\\n{2}{0}={1!r}'.format(key, value, indent))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Dumps representation to string. """
    from ..tools.uirepr import add_to_imports
    result = self.__repr_header__()
    indent = ' '.join('' for i in xrange(result.find('(')+1))
    add_to_imports(self, imports)

    # prints atoms.
    for o in self.atoms: 
      dummy = repr(o)
      dummy = dummy[dummy.find('(')+1:dummy.rfind(')')].rstrip().lstrip()
      result += '\\\n{0}.add_atom({1})'.format(indent, dummy)

    for item in self: 
      if item.__class__ is BaseKeyword:
        args = [repr(item.keyword)]
        if getattr(item, 'raw', None) is not None: args.append(repr(item.raw))
        result += '\\\n{0}.append({1})\n'.format(indent, ', '.join(args)) 
      else:
        add_to_imports(item, imports)
        item = repr(item).rstrip().lstrip()
        item = item.replace('\n', indent)
        result += '\\\n{0}.append({1})\n'.format(indent, item) 

    result = result.rstrip()
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())
    return {name: result.rstrip()}

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

  def read_input(self, tree, owner=None, **kwargs):
    """ Parses an input tree. """
    from ..tools.input import Tree
    from . import registered
    self[:] = []
    self.raw = tree.prefix
    do_breaksym = False
    has_breakkeep = False
    for key, value in tree:
      # parses sub-block.
      if isinstance(value, Tree): continue
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
      else: newobject = BaseKeyword(keyword=key)
      if len(value) > 0: 
        try: newobject.raw = getattr(value, 'raw', value)
        except: pass
      self.append(newobject)
      do_breaksym = False
      has_breakkeep = False

  def append(self, keyword, raw=None, breaksym=None):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    if breaksym is not None:
      self.append(GeomKeyword(keyword=keyword, raw=raw, breaksym=None))
    else: super(Molecule, self).append(keyword, raw)
    return self

  def eval(self):
    """ Evaluates current structure. 
    
        Runs crystal to evaluate current structure.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. error import GrepError
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

        try: return Extract().input_structure
        except GrepError:
          message = self.print_input() + '\n\n'
          with Extract().__stdout__() as file: message += file.read()
          raise ValueError(message + '\n\nInput tructure is likely incorrect\n')
    finally:
      try: rmtree(tmpdir)
      except: pass

  @property
  def crystal_output(self):
    """ Crystal output for the geometry alone. """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. import crystal_program as program
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
        with open('crystal.out', 'r') as file: return file.read()
    finally:
      try: rmtree(tmpdir)
      except: pass

  @property
  def symmetry_operators(self):
    """ Symmetry operators, as determined by CRYSTAL.
    
        Runs crystal to evaluate the symmetry operators on the current
        structure.
    """
    from copy import deepcopy
    from tempfile import mkdtemp
    from shutil import rmtree
    from ..misc import Changedir
    from ..process import ProgramProcess as Process
    from .. error import GrepError
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

        try: return Extract().symmetry_operators
        except GrepError:
          message = self.print_input() + '\n\n'
          with Extract().__stdout__() as file: message += file.read()
          raise ValueError(message + '\n\nInput tructure is likely incorrect\n')
    finally:
      try: rmtree(tmpdir)
      except: pass

  def copy(self):
    """ Returns deep copy of self. """
    from copy import deepcopy
    return deepcopy(self)

  def print_input(self, **kwargs):
    """ Prints as CRYSTAL output. """
    from .input import print_input
    return print_input(self.output_map(**kwargs))
