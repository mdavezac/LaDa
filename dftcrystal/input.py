__docformat__ = "restructuredtext en"
__all__ = ['GeomKeyword', 'AttrBlock']
from ..tools.input import AttrBlock as AttrBlockBase, BaseKeyword

class GeomKeyword(BaseKeyword):
  """ Keyword with breaksymm attribute. """
  def __init__(self, keyword=None, raw=None, **kwargs):
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
    self.breaksym = kwargs.get('breaksym', not kwargs.get('keepsym', True))
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
  def __repr__(self): 
    """ Dumps representation to string. """
    args = []
    if 'keyword' in self.__dict__:
      args.append("keyword={0.keyword!r}".format(self))
    if 'raw' in self.__dict__: args.append("raw={0.raw!r}".format(self))
    if self.breaksym == True: args.append("breaksym=True")
    return "{0.__class__.__name__}(".format(self) + ', '.join(args) + ')'
  def output_map(self, **kwargs):
    """ Print input to crystal. """
    from collections import OrderedDict
    # starts block
    result = OrderedDict(keepsymm = True) if self.keepsym                      \
             else OrderedDict(breaksym= True)
    # prints raw input, if present.
    raw = getattr(self, 'raw', None)
    if raw is not None:
      raw = raw.rstrip().lstrip()
      if len(raw) == 0: raw = None
    result[self.keywpord] = raw
    return result

class AttrBlock(AttrBlockBase):
  """ Defines block input to CRYSTAL. 
  
      A block is a any set of input which starts with a keyword and ends with
      an END. It can contain other sub-keywords.
      This particular flavor creates attributes from the  inner input keywords.
      This supposes that each keyword is only ever inputed once. 

      It can contain subitems.
  """
  def __init__(self, keyword=None, raw=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    super(AttrBlock, self).__init__(keyword=keyword, raw=raw)
    
  def output_map(self, **kwargs):
    """ CRYSTAL input.
    
        The input is composed as followed:
        
          - the keyword of this block, if it exists.
          - the 'raw' input, if is exists.
          - each input keyword to CRYSTAL contained by the block.
    """
    raw = getattr(self, 'raw', None)
    result = super(AttrBlock, self).output_map(**kwargs)
    if result is None and raw is None: return None
    if raw is not None: result.prefix = raw
    return result

  def read_input(self, tree, owner=None):
    """ Parses an input tree. """
    from ..error import internal
    from ..tools.input import Tree
    from . import registered
    if hasattr(self, 'raw') and len(getattr(tree, 'prefix', '')) > 0: 
      try: self.raw = tree.prefix
      except: pass
    do_breaksym = False
    has_breakkeep = False
    for key, value in tree:
      # parses sub-block.
      # check for special keywords.
      if key.lower() == 'keepsymm':
        do_breaksym, has_breakkeep = False, True
        continue
      if key.lower() == 'breaksym':
        do_breaksym, has_breakkeep = True, True
        continue
      if key.lower() in self._input: 
        newobject = self._input[key.lower()] 
        if hasattr(newobject, 'read_input'):
          newobject.read_input(value, owner=self)
        elif hasattr(newobject, 'raw'): newobject.raw = value
        elif hasattr(newobject, '__set__'): newobject.__set__(self, value)
        else:
          raise internal( "LaDa doesn't understand how to read input to {0}"   \
                          .format(key.lower()) )
        if has_breakkeep and hasattr(newobject, 'breaksym'):
          newobject.breaksym = do_breaksym
        do_breaksym, has_breakkeep = False, False
        continue
      if isinstance(value, Tree):
        newobject = registered.get(key.lower(), AttrBlock)()
        newobject.read_input(value, owner=self)
      # creates new object.
      elif key.lower() in registered:
        newobject = registered[key.lower()]()
        if len(value) > 0:
          try: newobject.raw = getattr(value, 'raw', value)
          except: pass
      elif len(value) == 0: newobject = True
      else: newobject = value
      self.add_keyword(key.lower(), newobject)
      has_breakkeep = False


def print_input(map):
  """ Prints output map to string. """
  from ..tools.input import Tree
  result = ''
  for key, item in map.iteritems():
    if isinstance(item, Tree):
      prefix = getattr(item, 'prefix', None)
      if hasattr(prefix, '__len__') and len(prefix) == 0: prefix = None
      result += key.upper()  + '\n'
      if prefix is not None: result += prefix.upper().rstrip() + '\n'
      dummy = print_input(item).rstrip().lstrip()
      if len(dummy) > 0: result += dummy.upper() + '\n'
      result += 'END ' + key.upper() + '\n'
    else: 
      if item in ['True', 'None']: 
        result += str(key).upper().rstrip().lstrip() + '\n'
      elif item is False: continue
      elif item is True or item is None:
        result += key.upper() + '\n'
      elif item is not None and len(item) > 0 and item != 'False':
        result += key.upper().rstrip().lstrip() + '\n'                         \
                  + item.upper().rstrip() + '\n'
  return result

  def _read_nested_group(self, tree, owner=None, **kwargs):
    """ Creates nested groups.

        It is not easy to guess the type of a yet to be created object.
        This method is used by derived classes to create default nested groups.
    """
    result = self.__class__()
    result.read_input(tree, owner=self, **kwargs)
    return result
