class Keyword(object):
  """ Defines keyword input to CRYSTAL. """
  def __init__(self, keyword=None, raw=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    super(Keyword, self).__init__()
    if keyword != None: 
      self.keyword = keyword
      """ CRYSTAL keyword. """
    if raw is not None:
      self.raw = raw
      """ Extra input to keyword. """

  def __repr__(self): 
    """ Dumps representation to string. """
    result = "{0.__class__.__name__}({0.keyword!r}".format(self)
    if 'raw' in self.__dict__: result += ", {0.raw!r}".format(self)
    return result + ')'
  
  def print_input(self, **kwargs):
    """ Print input to crystal. """
    # starts block
    result = '{0}\n'.format(self.keyword.uppper())

    # prints raw input, if present.
    raw = getattr(self, 'raw', None)
    if raw is not None:
      raw = self.raw.rstrip().lstrip()
      if len(raw) is not None: result += raw
    if result[-1] != '\n': result += '\n'
    return result

class GeomKeyword(Keyword):
  """ Adds breaksymm to :py:class:`~lada.dftcrystal.input.Keyword`. """
  def __init__(self, keyword=None, raw=None, *kwargs):
    """ Creates a geometry keyword. 
    
        :param str keyword: 
          keyword identifying the block.
        :param bool keepsym: 
          Whether to keep symmetries or not. Defaults to True.
        :param bool breaksym:
          Whether to break symmetries or not. Defaults to False.
          Only one of breaksymm needs be specified. If both are, then they
          should be consistent.
    """
    super(GeomKeyword, self).__init__(keyword=keyword, raw=raw)
    if 'keepsym' in kwargs and 'breaksym' in kwargs:
      if kwargs['keepsym'] == kwargs['breaksym']:
        raise ValueError('keepsym and breaksym both specified and equal.')
    self.breaksym = kwargs.get('breaksym', not kwargs.get('keepsym', False))
    """ Whether or not to break symmetries.
    
        Defaults to False. If symmetries are not broken, then all equivalent
        atoms are removed.
    """
    kwargs.pop('breaksym', None)
    kwargs.pop('keepsym', None)
  @property
  def keepsym(self):
    """ Not an alias for breaksym. """
    return not self.breaksym
  @keepsym.setter
  def keepsym(self, value):
    self.breaksym = not value

class Block(Keyword, list):
  """ Defines block input to CRYSTAL. 
  
      A block is a any set of input which starts with a keyword and ends with
      an END. It can contain other sub-keywords.

      It can contain subitems.
  """
  _blocks = [ 'crystal', 'slab', 'polymer', 'helix',  'molecule', 'freqcalc',
              'optgeom', 'anharm']
  """ Names of possible subblocks. 
  
      Used to parse input.
  """
  def __init__(self, keyword=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    Keyword.__init__(self, keyword)
    list.__init__(self)

  def append(self, value):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    list.append(self, value)
    return self

  def __repr__(self, indent=''): 
    """ Dumps representation to string. """
    from inspect import getargspec
    result = indent + Keyword.__repr__(self)
    if len(self) > 0:
      indent += '    '
      for item in self: 
        hasindent = getargspec(item.__repr__).keywords
        hasindent = False if hasindent is None else 'indent' in hasindent
        result += '\\\n{0},append{1!r}'.format(indent, item) if not hasindent  \
                  else '\\\n{0}{1}'.format(indent, item.__repr__(indent)) 
    if len(self) > 0: result = result[:result.rfind(')')+1]
    return result
  
  def print_input(self, **kwargs):
    """ Print input to crystal. """
    # starts block
    result = Keyword.print_input(self, **kwargs)

    # adds subtitems if present.
    if len(self._operators) > 0:
      for op in self._operators: 
        dummy = op.print_input(**kwargs) 
        if dummy is None: continue
        dummy = dummy.lstrip()
        if dummy[-1] != '\n': dummy += '\n'
        result += dummy

    # ends block
    result += 'END {0}'.format(self.keyword.upper())
    return result

  def find_last(self, key):
    """ Finds last of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    return self[items.rfind(key)]
  def find_first(self, key):
    """ Finds first of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    return self[items.find(key)]

  def pop_last(self, key):
    """ Pops last of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    return self.pop(items.rfind(key))
  def pop_first(self, key):
    """ Pops first of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    return self.pop(items.rfind(key))

  def replace_last(self, key, value):
    """ Replaces last of kind. """
    from ..error import KeyError
    if not (isinstance(value, Keyword) or isinstance(value, Block)): 
      raise ValueError('Input is not a crystal Keyword or Block')
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    self[items.rfind(key)] = value
  def replace_first(self, key, value):
    """ Replaces first of kind. """
    from ..error import KeyError
    if not (isinstance(value, Keyword) or isinstance(value, Block)): 
      raise ValueError('Input is not a crystal Keyword or Block')
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    self[items.find(key)] = value

class Input(list):
  __slots__ = ['raw']
  def __init__(self):
    super(Input, self).__init__()
  def descend(self, *args):
    if len(args) == 0: return self
    name, args = args[0], args[1:]
    for key, value in self: 
      if name == key:
        return value.descend(*args) if hasattr(value, 'descend') else value
    self.append((name, Input()))
    return self[-1][1].descend(*args)
  def __getitem__(self, name):
    if isinstance(name, str): 
      for key, value in self: 
        if name == key: return value
    return super(Input, self).__getitem__(name)
  def keys(self): return [u[0] for u in self]
  def __contains__(self, value): return value in self.keys()
    
  
def parse_input(path):
  """ Reads crystal input. """
  from re import compile
  from copy import copy
  if isinstance(path, str): 
    with open(path) as file: return parse_input(file)

  starters = ['CRYSTAL', 'SLAB', 'POLYMER', 'HELIX', 'MOLECULE']
  title = ''
  for i, line in enumerate(path):
    if line.split()[0] in starters: break
    title = line
  keyword_re = compile('^[A-Z](?!\s)')

  blocks = copy(starters)
  blocks.extend(['OPTGEOM', 'FREQCALC', 'ANHARM', 'DFT'])
  
  # reading linearly, 
  title = title[:-1].rstrip().lstrip()
  nesting = [title, line.split()[0]]
  results = Input()
  keyword = line.split()[0]
  raw = ''
  # reads crystal input.
  for line in path:
    # found end, pop nesting.
    if line.split()[0][:3] == 'END': 
      current = nesting.pop(-1)
      if current in starters: 
        nesting.append('BASISSET')
      if len(nesting) == 0: break
    # special case of INPUT keyword
    elif line.split()[0] == 'INPUT': raw += line
    # found keyword
    elif keyword_re.match(line) is not None:
      newkeyword = line.split()[0]
      current = results.descend(*nesting)
      # first of subblock
      if keyword == nesting[-1]: current.raw = raw
      if newkeyword in blocks                                                  \
         and not (newkeyword == 'SLAB' and nesting[-1] == 'CRYSTAL'): 
        nesting.append(newkeyword)
      # normal keyword
      else: current.append((newkeyword, raw)) 
      raw = ''
      keyword = newkeyword
    # found raw string
    else: raw += line
  return results
