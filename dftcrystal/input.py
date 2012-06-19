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

class Block(Keyword):
  """ Defines block input to CRYSTAL. 
  
      A block is a any set of input which starts with a keyword and ends with
      an END. It can contain other sub-keywords.

      It can contain subitems.
  """
  def __init__(self, keyword):
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
