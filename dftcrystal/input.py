__docformat__ = "restructuredtext en"
__all__ = ['Keyword', 'GeomKeyword', 'ListBlock', 'ValueKeyword']
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
    args = []
    if 'keyword' in self.__dict__:
      args.append("keyword={0.keyword!r}".format(self))
    if 'raw' in self.__dict__: args.append("raw={0.raw!r}".format(self))
    return "{0.__class__.__name__}(".format(self) + ', '.join(args) + ')'
  
  def print_input(self, **kwargs):
    """ Print input to crystal. """
    # starts block
    result = '{0}\n'.format(self.keyword.upper())

    # prints raw input, if present.
    raw = getattr(self, 'raw', None)
    if raw is not None:
      raw = self.raw.rstrip().lstrip()
      if len(raw) is not None: result += raw
    if result[-1] != '\n': result += '\n'
    return result

class ValueKeyword(Keyword):
  """ Keywords which expect a value of some kind. 
  
  
      LaDa will attempt to guess that value when set using
      :py:attr:`~ValueKeyword.raw`. It will use it as is when set through
      :py:meth:`~ValueKeyword.__set__`. Finally, it will not print if the value
      is None.
  """
  def __init__(self, keyword, value=None):
    """ Initializes a keyword with a value. """
    super(ValueKeyword, self).__init__(keyword=keyword)
    self.value = value
    """ The value to print to input. 

        If None, then this keyword will not appear.
    """
  def print_input(self, **kwargs):
    """ Prints keyword to string. """
    if self.value == None: return None
    return '{0}\n{1}\n'.format(self.keyword.upper(), self.raw)
  @property
  def raw(self):
    """ Returns raw value for CRYSTAL input. """
    return str(self.value).upper() if not hasattr(self.value, '__iter__')      \
           else ' '.join(str(u).upper() for u in self.value)
  @raw.setter
  def raw(self, value):
    """ Guesses value from raw input. """
    from ..error import ValueError
    if not isinstance(value, str):
      raise ValueError('Expected a string as input to raw.')
    # split along line and remove empty lines
    lines = value.split('\n')
    while len(lines[-1].rstrip().lstrip()) == 0: lines.pop(-1)
    # if only one line left, than split into a list and guess type of each
    # element.
    if len(lines) == 0:
      raise ValueError('Expected *non-empty* string as input.\n')
    elif len(lines) == 1:
      lines = lines[0]
      n = []
      for u in lines.split():
        try: v = int(u)
        except:
          try: v = float(u)
          except: v = u
        n.append(v)
      # if only one element use that rather than list
      if len(n) == 1: n = n[0]
    # if multiple line, keep as such
    else: n = value
    self.value = n
  def __get__(self, instance, owner=None):
    """ Returns value as is. """
    return self.value
  def __set__(self, instance, value):
    """ Assigns value as is. """
    self.value = value
  

class TypedKeyword(Keyword):
  """ Keywords which expect a value of a given type.
  
  
      Two types are allowed: 
      
        - None, in which case the keyword won't appear in the input
        - The explicitely given type, eg int.
        - A list of classes. Yes, a list, not a tuple or any other sequence.

            - only one item: the value can have any length, but the
              type of each element must conform. For example ``[int]`` will map
              "5 6 7" to ``[5, 6, 7]``. 
            - more than one item: the value is a list n items, where n is the
              size of the type. Each element in the value must conform to the
              respective element in the type. For example, ``[int, str]`` will
              map "5 two" to ``[5, "two"]. It will fail when mapping "5 6 7"
              (three elements) or "two 5" (wrong types).

              .. warning:: 

                 Types are checked if the value is set as a whole, not when a
                 single item is set.
     
      The value given will be cast to the given type, unless None, in which
      case the keyword will not be printed to the CRYSTAL input file.
  """
  def __init__(self, keyword=None, type=int, value=None):
    """ Initializes a keyword with a value. """
    from ..error import ValueError
    super(TypedKeyword, self).__init__(keyword=keyword)
    if isinstance(type, list) and len(type) == 0:
      raise ValueError('type must be class or a non-empty list of classes')

    self.type = type
    """ Type to which the value should be cast if the value is not None. """
    self.value = value
    """ The value to print to input. 

        If None, then this keyword will not appear.
    """
  @property
  def value(self): 
    """ The value to print to input. 

        If None, then this keyword will not appear.
        Otherwise, it is cast to the type given when initializing this instance
        of :py:class:`~lada.dftcrystal.input.TypedKeyword`
    """
    return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return
    if type(self.type) is list:
      if not hasattr(value, '__iter__'): 
        raise ValueError('Expected a sequence on input.')
      if len(self.type) == 1: 
        _type = self.type[0]
        self._value = [_type(u) for u in value]
        if len(self._value) == 0: self._value = None
      else: 
        if len(value) != len(self.type):
          raise ValueError( 'Expected a sequence of the following type: {0}'   \
                            .format(self.type) )
        self._value = [t(v) for t, v in zip(self.type, value)]
    else: self._value = self.type(value)
  def print_input(self, **kwargs):
    """ Prints keyword to string. """
    if self._value == None: return None
    return '{0}\n{1}\n'.format(self.keyword.upper(), self.raw)
  @property
  def raw(self):
    """ Returns raw value for CRYSTAL input. """
    if type(self.type) is list:
      return ' '.join(str(v).upper() for v in self.value)
    return str(self._value).upper()
  @raw.setter
  def raw(self, value):
    """ Guesses value from raw input. """
    if type(self.type) is list: value = value.split()
    self.value = value
  def __get__(self, instance, owner=None):
    """ Returns value as is. """
    return self.value
  def __set__(self, instance, value):
    """ Assigns value as is. """
    self.value = value

class BoolKeyword(Keyword):
  """ Boolean keyword.

      If True, the keyword is present.
      If False, it is not.
      This class uses the get/set mechanism to set whether the keyword should
      appear or not. It is meant to be used in conjunction with other linked keywords.
      Otherwise, it is simpler to use :py:meth:`self.add_keyword('something')
      <lada.dftcrystal.input.AttrBlock.add_keyword>` directly.

      .. note:: 

        This class should be derived from. It is not complete since it lacks a
        keyword attribute.
  """
  def __init__(self, value=False, keyword=None):
    """ Initializes FullOptG keyword. """
    super(BoolKeyword, self).__init__(keyword=keyword)
    self.value = value
  @property
  def value(self): return self._value
  @value.setter
  def value(self, value): self._value = (value == True)
  def __set__(self, instance, value):
    """ Sets the keyword to appear or not. """
    self.value = value
  def __get__(self, instance, owner=None):
    """ True if the keyword is to appear. """
    return self.value
  def print_input(self, **kwargs):
    return self.keyword.upper() + '\n' if self.value else None

class GeomKeyword(Keyword):
  """ Adds breaksymm to :py:class:`~lada.dftcrystal.input.Keyword`. """
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
  def print_input(self, **kwargs):
    """ Print input to crystal. """
    # starts block
    result = 'KEEPSYMM\n' if self.keepsym else 'BREAKSYM\n'
    result += '{0}\n'.format(self.keyword.upper())

    # prints raw input, if present.
    raw = getattr(self, 'raw', None)
    if raw is not None:
      raw = self.raw.rstrip().lstrip()
      if len(raw) is not None: result += raw
    if result[-1] != '\n': result += '\n'
    return result

class ListBlock(Keyword, list):
  """ Defines block input to CRYSTAL. 
  
      A block is a any set of input which starts with a keyword and ends with
      an END. It can contain other sub-keywords.
      This particular flavor appends inner input keywords to a list. It is used
      within the geometry block of CRYSTAL, since that particular block is list
      of commands which are executed one after the other. 

      It can contain subitems.
  """
  def __init__(self, keyword=None, raw=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    Keyword.__init__(self, keyword=keyword, raw=raw)
    list.__init__(self)

  def append(self, keyword, raw=None):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    from ..error import ValueError
    if isinstance(keyword, str):
      keyword = Keyword(keyword=keyword, raw=raw)
    elif not isinstance(keyword, Keyword):
      raise ValueError('Wrong argument to append.')
    list.append(self, keyword)
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
    for op in self: 
      dummy = op.print_input(**kwargs) 
      if dummy is None: continue
      dummy = dummy.lstrip()
      if dummy[-1] != '\n': dummy += '\n'
      result += dummy

    # ends block
    result += 'END {0}\n'.format(self.keyword.upper())
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
    if not isinstance(value, Keyword):
      raise ValueError('Input is not a crystal Keyword or Block')
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    self[items.rfind(key)] = value
  def replace_first(self, key, value):
    """ Replaces first of kind. """
    from ..error import KeyError
    if not isinstance(value, Keyword):
      raise ValueError('Input is not a crystal Keyword or Block')
    items = [u.keyword for u in self] 
    if key not in items: raise KeyError('Could not find {0}.'.format(key))
    self[items.find(key)] = value

class AttrBlock(Keyword):
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
    # first add this to avoid infinite recursion from redifining __setattr__
    self.__dict__['_crysinput'] = {}

    # then call base constructor
    super(AttrBlock, self).__init__(keyword=keyword, raw=raw)

    # now so we get doctstrings right.
    self._crysinput = {}
    """ Dictionary of crystal inputs. """
    

  def __getattr__(self, name):
    """ passes through the input keywords in :py:attr:`_crysinput`. """
    from ..error import AttributeError
    if name not in self._crysinput: 
      raise AttributeError('Unknwon attribute {0}.'.format(name))
    result = self._crysinput[name]
    return result.__get__(self) if hasattr(result, '__get__') else result
  def __setattr__(self, name, value):
    """ passes through the input keywords in :py:attr:`_crysinput`. 
    
        If the input value is derived from
        :py:class:`~lada.dftcrystal.input.Keyword`, then it is added to
        :py:attr:`_crysinput`. Otherwise, super is called.
    """
    if isinstance(value, Keyword): 
      self._crysinput[name] = value
    elif name in self._crysinput:
      result = self._crysinput[name]
      if hasattr(result, '__set__'): result.__set__(self, value)
      else: self._crysinput[name] = value
    else: super(AttrBlock, self).__setattr__(name, value)
  def __delattr__(self, name):
    """ passes through the input keywords in :py:attr:`_crysinput`.  """
    if name in self._crysinput: del self._crysinput[name]
    else: super(AttrBlock, self).__delattr__(name)

  def add_keyword(self, name, value=None):
    """ Adds/Sets input keyword. """
    # if known keyword, then go through setattr mechanism.
    # this makes sure we recognize the type of value and the already registered
    # keyword.
    if name in self._crysinput:  setattr(self, name, value)
    # if value is None, then transform it to True. 
    # This is a keyword which is either there or not there, like EXTPRT.
    elif value is None: self._crysinput[name] = True
    # boolean case
    elif value is True or value is False:
      self._crysinput[name] = value
    # if a string, tries to guess what it is.
    elif isinstance(value, str):
      # split along line and remove empty lines
      lines = value.split('\n')
      while len(lines[-1].rstrip().lstrip()) == 0: lines.pop(-1)
      # if only one line left, than split into a list and guess type of each
      # element.
      if len(lines) == 1:
        lines = lines[0]
        n = []
        for u in lines.split():
          try: v = int(u)
          except:
            try: v = float(u)
            except: v = u
          n.append(v)
        # if only one element use that rather than list
        if len(n) == 1: n = n[0]
      # if multiple line, keep as such
      else: n = value
      self._crysinput[name] = n
    # otherwise, just set the keyword.
    else: self._crysinput[name] = value
    # return self to allow chaining calls.
    return self


  def print_input(self, **kwargs):
    """ CRYSTAL input.
    
        The input is composed as followed:
        
          - the keyword of this block, if it exists.
          - the 'raw' input, if is exists.
          - each input keyword to CRYSTAL contained by the block.
    """
    result = getattr(self, 'raw', '')
    for key, value in self._crysinput.iteritems():
      if value is None: continue
      elif isinstance(value, bool):
        if value == True: result += key.upper()
      elif hasattr(value, 'print_input'):
        dummy = value.print_input(**kwargs)
        if dummy is not None: result += dummy
      elif getattr(value, 'raw', None) is not None:
        result += getattr(result, 'keyword', key).upper() + '\n'               \
                  + str(value.raw)
      elif hasattr(value, '__iter__'):
        result += getattr(result, 'keyword', key).upper() + '\n'               \
                  + ' '.join(str(u) for u in value)
      else:
        result += getattr(result, 'keyword', key).upper() + '\n' + str(value)
      result = result.rstrip()
      if len(result) > 0 and result[-1] != '\n': result += '\n'
    if len(result) == 0: return None # nothing printed, return None.
    if result[-1] != '\n': result += '\n'
    if getattr(self, 'keyword', '') != '': 
      return '{0}\n{1}END {0}\n'.format(self.keyword.upper(), result)
    return result

  def __repr__(self, indent=''):
    """ Representation of this instance. """
    from inspect import getargspec
    result = indent + '{0.__class__.__name__}()'.format(self)
    indent += '    '
    for key, value in self._crysinput.iteritems():
      if hasattr(getattr(value, '__repr__', None), '__call__'):
        try:
          hasindent = getargspec(value.__repr__).args
        except: hasindent = False
        else: hasindent = False if hasindent is None else 'indent' in hasindent
      else: hasindent = False
      if hasindent: 
        result += '\\\n{0}.add_keyword({1!r}, {2})'                            \
                  .format(indent, key, value.__repr__(indent+'        ')) 
      else: 
        result += '\\\n{0}.add_keyword({1!r}, {2!r})'                          \
                  .format(indent, key, value) 
    return result

  def read_input(self, tree):
    """ Parses an input tree. """
    from parse import InputTree
    from . import registered
    if hasattr(self, 'raw') and len(tree.raw) > 0: 
      try: self.raw = tree.raw
      except: pass
    for key, value in tree:
      # parses sub-block.
      if key.lower() in self._crysinput: 
        newobject = self._crysinput[key.lower()] 
        if hasattr(newobject, 'read_input'):
          newobject.read_input(value)
          continue
        elif hasattr(newobject, '__set__'):
          newobject.__set__(self, value)
        elif hasattr(newobject, 'raw'):
          newobject.raw = value
          continue
      if isinstance(value, InputTree):
        newobject = registered.get(key.lower(), AttrBlock)()
        newobject.read_input(value)
      # creates new object.
      elif key.lower() in registered:
        newobject = registered[key.lower()]()
        if len(value) > 0:
          try: newobject.raw = getattr(value, 'raw', value)
          except: pass
      elif len(value) == 0: newobject = True
      else: newobject = value
      self.add_keyword(key.lower(), newobject)


class Choice(Keyword):
  """ Keyword value must be chosen from a given set. """
  def __init__(self, values=None, value=None, keyword=None):
    """ Creates keyword which must be chosen from a given set. """ 
    super(Choice, self).__init__(keyword=keyword)
    if values is not None:
      self.values = list(values)
      """ Set of values from which to choose keyword. """
    self.value = value
    """ Current value. """
  @property
  def value(self):
    """ Current value of the keyword. """
    return self._value
  @value.setter
  def value(self, value):
    from ..error import ValueError
    if value is None: self._value = None; return
    if hasattr(value, 'rstrip'): value = value.rstrip().lstrip()
    if hasattr(value, 'lower'): value = value.lower()
    for v in self.values:
      try: o = v.__class__(value)
      except: pass
      if (hasattr(o, 'lower') and o.lower() == v.lower()) or o == v: 
        self._value = o
        return
    raise ValueError( 'Input({0!r}) should be one of the following: {1}.'      \
                      .format(value, self.values) )
  def __get__(self, instance, owner=None):
    """ Function called by :py:class:`AttrBlock`. """
    return self._value
  def __set__(self, instance, value):
    """ Function called by :py:class:`AttrBlock`. """
    self.value = value

  @property
  def raw(self): return str(self.value)
  @raw.setter
  def raw(self, value): self.value = value
    
  def print_input(self, **kwargs):
    """ Prints input to string. """
    if self._value is None: return None
    result = getattr(self, 'keyword', '').upper()
    if len(result) > 0 and result[-1] != '\n': result += '\n'
    result += str(self.value).upper()
    return result + '\n'
