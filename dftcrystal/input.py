__docformat__ = "restructuredtext en"
__all__ = [ 'Keyword', 'GeomKeyword', 'ListBlock', 'ValueKeyword',
            'TypedKeyword', 'VariableListKeyword', 'QuantityKeyword' ]
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
  
      Instances of this class make it easy to declare and use CRYSTAL keywords
      which define a single value::

        functional.keyword = ValueKeyword(value=5)

      The above will print to the Crystal input as follows:

        | KEYWORD
        | 5

      The keyword can be then be given any value::

        functional.keyword = abc
      
      would result in:

        | KEYWORD
        | ABC

      In practice, the value is printed by first transforming it to a string
      via str_ and putting it in upper case.

      If given a string, LaDa will attempt to guess the type of the value,
      depending on whether it is an ``int``, ``float``, ``str``, or a list of
      the same, in that order. A list is created if and only if the string is
      does not contian multiple lines. If the string contains more than one
      line, it will always be kept as a string, so that complex pre-formatted
      input is not messed with. :py:meth:`~ValueKeyword.__set__`. Finally, the
      keyword will not appear in the input if its value is None. 

      .. _str: http://docs.python.org/library/functions.html#str
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
    return '{0}\n{1}\n'.format(self.keyword.upper(), self.raw.rstrip())
  @property
  def raw(self):
    """ Returns raw value for CRYSTAL input. """
    if self.value == None: return '' # otherwise, fails to find attribute.
    return str(self.value).upper() if not hasattr(self.value, '__iter__')      \
           else ' '.join(str(u).upper() for u in self.value)
  @raw.setter
  def raw(self, value):
    """ Guesses value from raw input. """
    from ..error import ValueError
    if not isinstance(value, str):
      raise ValueError( 'Expected a string as input to {0}.raw.'               \
                        .format(self.keyword) )
    # split along line and remove empty lines
    lines = value.split('\n')
    while len(lines[-1].rstrip().lstrip()) == 0: lines.pop(-1)
    # if only one line left, than split into a list and guess type of each
    # element.
    if len(lines) == 0:
      raise ValueError( 'Expected *non-empty* string as '                      \
                        'input to {0.keyword}.\n'.format(self) )
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
  def __repr__(self):
    """ Dumps a representation of self. """
    from inspect import getargspec
    args = []
    if 'keyword' not in self.__class__.__dict__ and 'keyword' in self.__dict__:
      args.append('keyword={0.keyword!r}'.format(self))
    argspec = getargspec(self.__init__)
    index = argspec.args.index('value')                                         \
            - len(argspec.args) + len(argspec.defaults)
    dovalue = True
    if index >= 0:
      default = argspec.defaults[index]
      if default is None and self.value is None: dovalue = False
      else:
        try: 
          if default.__class__(self.value) == default: dovalue = False
        except: pass
    if dovalue: args.append('value={0.value!r}'.format(self))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))
  

class TypedKeyword(ValueKeyword):
  """ Keywords which expect a value of a given type.
  
      This specializes :py:class:`ValueKeyword` to accept only input of a given
      type.
  
      Two types are allowed: 
      
        - A simple type, eg int.
        - A list of types. Yes, a list, not a tuple or any other sequence.

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

      In practice, we can create the following mapping from python to the
      CRYSTAL input::

         functional.keyword = TypedValue([int, float, int])
         functional.keyword = [5, 6, 3]

      will yield:

         | KEYWORD
         | 5 6.0 3

      The following would throw an exception:

        >>> functional.keyword = ['a', 5.3, 5]
        ValueError: ...
        >>> functional.keyword = [0, 2.0]
        lada.error.ValueError: Expected a sequence of the following type: [int, float, int]

      The last case failed because only two values are given. For a list of
      variable length (but only one type of lement), use::

        functional.keyword = TypedKeyword(type=[int])
        functional.keyword = [0, 1, 2]
        functional.keyword = [0, 1, 2, 3]

      The last two lines will not raise an exception. The CRYSTAL input would
      look like this:

        | KEYWORD
        | 0 1 2 3 

      Note that the first value is *not* the length of the list! If you need
      that behavior, see :py:class:`VariableListKeyword`. 

      .. note::

        Formally, the type need not be set in
        :py:method:`~TypedKeyword.__init__`. If it is left out (in which case
        it defaults to None), it is expected that it will be set later by the
        user, prior to use. 
        This is mainly to allow derived classes to define
        :py:attr:`~TypedKeyword.type` as a class attribute, rather than an
        instance attribute.
  """
  def __init__(self, keyword=None, type=None, value=None):
    """ Initializes a keyword with a value. """
    from ..error import ValueError
    super(TypedKeyword, self).__init__(keyword=keyword)
    if isinstance(type, list) and len(type) == 0:
      raise ValueError('type must be class or a non-empty list of classes')

    if type is not None:
      self.type = type
      """ Type to which the value should be cast if the value is not None. """
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
        raise ValueError( '{0} expected a sequence on input, got {1!r}.'       \
                          .format(self.keyword, value) ) 
      if len(self.type) == 1: 
        _type = self.type[0]
        self._value = [_type(u) for u in value]
        if len(self._value) == 0: self._value = None
      else: 
        if len(value) != len(self.type):
          raise ValueError( '{0.keyword} expected a sequence of the '          \
                            'following type: {0.type}'.format(self) )
        self._value = [t(v) for t, v in zip(self.type, value)]
    else: self._value = self.type(value)
  @property
  def raw(self):
    """ Returns raw value for CRYSTAL input. """
    if self._value == None: return '' # otherwise, fails to find attribute.
    if type(self.type) is list:
      return ' '.join(str(v).upper() for v in self.value)
    return str(self._value).upper()
  @raw.setter
  def raw(self, value):
    """ Guesses value from raw input. """
    if type(self.type) is list: value = value.split()
    self.value = value
  def __repr__(self):
    """ Dumps a representation of self. """
    from inspect import getargspec
    args = []
    if 'keyword' not in self.__class__.__dict__ and 'keyword' in self.__dict__:
      args.append('keyword={0.keyword!r}'.format(self))
    if 'type' not in self.__class__.__dict__:
      if isinstance(self.type, list):
        args.append( 'type=[{0}]'                                              \
                     .format(', '.join(u.__name__ for u in self.type)) )
      else: args.append( 'type={0.type.__name__}'.format(self))
    argspec = getargspec(self.__init__)
    index = argspec.args.index('value') - 1
    dovalue = True
    if len(argspec.defaults) > index:
      default = argspec.defaults[index]
      if default is None and self.value is None: dovalue = False
      else:
        try: 
          if default.__class__(self.value) == default: dovalue = False
        except: pass
    if dovalue: args.append('value={0.value!r}'.format(self))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

class VariableListKeyword(TypedKeyword):
  """ Keywords which expect a variable-length list value of a given type.
  
      Expects a  list of values of a given type.
      The length can be anything. However, unlike its base class
      :py:class:`TypedKeyword`, the fortran input
      :py:attr:`~VariableListKeyword.raw` should consist first of an integer
      giving the size of the list of values.

      Furthermore, the type of the element of the list should be given on
      initialization. Eg, ``type=int``, rather than, say, ``type=[int]``
  """
  def __init__(self, keyword=None, type=None, value=None):
    """ Initializes a keyword with a value. """
    super(VariableListKeyword, self).__init__( keyword=keyword, type=type,     \
                                               value=value )

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
    if not hasattr(value, '__iter__'): 
      raise ValueError( '{0.keyword} expected a sequence on input.'            \
                        .format(self) )
    self._value = [self.type(u) for u in value]
    if len(self._value) == 0: self._value = None
  @property
  def raw(self):
    """ Returns raw value for CRYSTAL input. """
    if self._value == None: return '' # otherwise, fails to find attribute.
    lstr = ' '.join(str(v).upper() for v in self.value) 
    return '{0}\n{1}'.format(len(self.value), lstr).upper()
  @raw.setter
  def raw(self, value):
    """ Guesses value from raw input. """
    value = value.split()
    self.value = value[1:int(value[0])+1]


class BoolKeyword(ValueKeyword):
  """ Boolean keyword.

      If True, the keyword is present.
      If False, it is not.
      This class uses the get/set mechanism to set whether the keyword should
      appear or not. It is meant to be used in conjunction with other linked keywords.
      Otherwise, it is simpler to use :py:meth:`self.add_keyword('something')
      <lada.dftcrystal.input.AttrBlock.add_keyword>` directly.
  """
  def __init__(self, keyword=None, value=False):
    """ Initializes FullOptG keyword. """
    super(BoolKeyword, self).__init__(keyword=keyword, value=value)
  @property
  def value(self): return self._value
  @value.setter
  def value(self, value): self._value = (value == True)
  def print_input(self, **kwargs):
    return self.keyword.upper() + '\n' if self.value else None
  def read_input(self, tree, owner=None):
    """ Is True if ever read. """
    self.value = True

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
      raise ValueError('Wrong argument to {0.keyword}.append().'.format(self))
    list.append(self, keyword)
    return self

  def insert(self, index, keyword, raw=None):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    from ..error import ValueError
    if isinstance(keyword, str):
      keyword = Keyword(keyword=keyword, raw=raw)
    elif not isinstance(keyword, Keyword):
      raise ValueError('Wrong argument to {0.keyword}.append().'.format(self))
    list.insert(self, index, keyword)
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
    result = super(ListBlock, self).print_input(**kwargs)

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
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
    return self[items.rfind(key)]
  def find_first(self, key):
    """ Finds first of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
    return self[items.find(key)]

  def pop_last(self, key):
    """ Pops last of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
    return self.pop(items.rfind(key))
  def pop_first(self, key):
    """ Pops first of kind. """
    from ..error import KeyError
    items = [u.keyword for u in self] 
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
    return self.pop(items.rfind(key))

  def replace_last(self, key, value):
    """ Replaces last of kind. """
    from ..error import KeyError
    if not isinstance(value, Keyword):
      raise ValueError( 'Expected crystal Keyword as input to '                \
                        '{0.keyword}.replace_last()'.format(self) )
    items = [u.keyword for u in self] 
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
    self[items.rfind(key)] = value
  def replace_first(self, key, value):
    """ Replaces first of kind. """
    from ..error import KeyError
    if not isinstance(value, Keyword):
      raise ValueError( 'Expected crystal Keyword as input to '                \
                        '{0.keyword}.replace_last()'.format(self) )
    items = [u.keyword for u in self] 
    if key not in items:
      raise KeyError('Could not find {0} in {1.keyword}.'.format(key, self))
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
      raise AttributeError('Unknown attribute {0}.'.format(name))
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
      if not hasattr(value, 'keyword'): self._crysinput[name].keyword = name
    elif name in self._crysinput:
      result = self._crysinput[name]
      if hasattr(result, '__set__'): result.__set__(self, value)
      else: self._crysinput[name] = value
    else: super(AttrBlock, self).__setattr__(name, value)
  def __delattr__(self, name):
    """ passes through the input keywords in :py:attr:`_crysinput`.  """
    if name in self._crysinput: del self._crysinput[name]
    else: super(AttrBlock, self).__delattr__(name)
  def __dir__(self):
    """ List of attributes and members. """
    return list( set(self._crysinput.iterkeys())                               \
                 | set(self.__dict__.iterkeys())                               \
                 | set(dir(self.__class__)) )

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
      try: 
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
      except:
        from sys import exc_info
        type, value, traceback = exc_info()
        message = 'ERROR in when printing keyword {0}'.format(key)
        if value is None: type.args = type.args, message
        else: value = value, message
        raise type, value, traceback
    result = result.rstrip()
    if getattr(self, 'keyword', '') != '': 
      return '{0}\n{1}\nEND {0}\n'.format(self.keyword.upper(), result)
    return result

  def __repr__(self, defaults=False, name=None):
    """ Representation of this instance. """
    from ..functools.uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def read_input(self, tree, owner=None):
    """ Parses an input tree. """
    from ..error import internal
    from .parse import InputTree
    from . import registered
    if hasattr(self, 'raw') and len(tree.raw) > 0: 
      try: self.raw = tree.raw
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
      if key.lower() in self._crysinput: 
        newobject = self._crysinput[key.lower()] 
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
      if isinstance(value, InputTree):
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

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    from ..functools.uirepr import template_ui_repr, add_to_imports

    results = template_ui_repr(self, imports, name, defaults, exclude)
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())

    for key, value in self._crysinput.iteritems():
      if exclude is not None and key in exclude: continue
      if hasattr(value, '__ui_repr__'): 
        default = None if defaults is None                                     \
                  else defaults._crysinput.get(key, None)
        newname = name + '.' + key
        partial = value.__ui_repr__(imports, newname, default)
        results.update(partial)
        if newname in results:              doinit = False
        elif default is None:               doinit = True
        else: doinit = type(value) is not type(default)
        if doinit:
          results[newname] = '{0.__class__.__name__}()'.format(value)
          add_to_imports(value, imports)
      elif isinstance(value, Keyword):
        value = getattr(self, key)
        string = repr(value)
        if defaults is not None and key in defaults._crysinput                 \
           and type(value) is type(getattr(defaults, key))                     \
           and string == repr(getattr(defaults, key)): continue
        key = '{0}.{1}'.format(name, key) 
        results[key] = string
        add_to_imports(value, imports)
      elif value is None:
        if defaults is not None and key in defaults._crysinput                 \
           and defaults._crysinput[key] is None: continue
        results['{0}.add_keyword({1!r})'.format(name, key)] = None
      else:
        if defaults is not None and key in defaults._crysinput                 \
           and type(value) is type(defaults._crysinput[key])                   \
           and repr(value) == repr(defaults._crysinput[key]): continue
        results['{0}.add_keyword({1!r}, {2!r})'.format(name, key, value)]      \
            = None
        add_to_imports(value, imports)
    
    return results
  
  def __getstate__(self):
    d = self.__dict__.copy()
    crysinput = d.pop('_crysinput')
    return d, crysinput
  def __setstate__(self, value):
    self.__dict__['_crysinput'] = value[1]
    self.__dict__.update(value[0])

class ChoiceKeyword(Keyword):
  """ Keyword value must be chosen from a given set. """
  def __init__(self, values=None, value=None, keyword=None):
    """ Creates keyword which must be chosen from a given set. """ 
    super(ChoiceKeyword, self).__init__(keyword=keyword)
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
    raise ValueError( '{0.keyword} accepts only one of the following: {1}'     \
                      .format(self, self.values) )
  def __get__(self, instance, owner=None):
    """ Function called by :py:class:`AttrBlock`. """
    return self._value
  def __set__(self, instance, value):
    """ Function called by :py:class:`AttrBlock`. """
    self.value = value

  @property
  def raw(self):
    if self._value == None: return '' # otherwise, fails to find attribute.
    return str(self.value)
  @raw.setter
  def raw(self, value): self.value = value
    
  def print_input(self, **kwargs):
    """ Prints input to string. """
    if self._value is None: return None
    result = getattr(self, 'keyword', '').upper()
    if len(result) > 0 and result[-1] != '\n': result += '\n'
    result += str(self.value).upper()
    return result + '\n'

class QuantityKeyword(ValueKeyword):
  """ Keyword with a value which is signed by a unit. """
  def __init__(self, units=None, shape=None, keyword=None, value=None):
    """ Creates the quantity itself. """

    super(QuantityKeyword, self).__init__(keyword=keyword)
    if units is not None:
      self.units = units
      """ UnitQuantity to which values should be scaled. """
    if shape is not None:
      self.shape = shape
      """ Shape of input/output arrays. """
    self.value = value

  @property
  def value(self): 
    """ Value to which this keyword is set. """
    return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return
    if hasattr(value, 'rescale'): value = value.rescale(self.units)
    else: value = value * self.units
    shape = getattr(self, 'shape', ())
    if value.shape != shape: value = value.reshape(shape)
    self._value = value
  @property
  def raw(self):
    """ Returns string for CRYSTAL input. """
    if self._value is None: return ''
    shape = getattr(self, 'shape', ())
    if len(shape) > 2: 
      from ..error import NotImplementedError
      raise NotImplementedError( 'LaDa does not know how to print n-d arrays ' \
                                 '(n>2) to CRYSTAL input.')
    if len(shape) == 0: return str(float(self.value))
    else:
      result = str(self.value.magnitude).replace('[', ' ').replace(']', ' ')
      result = result.split('\n')
      return '\n'.join(u.rstrip().lstrip() for u in result)
  @raw.setter
  def raw(self, value):
    """ Creates value from CRYSTAL input. """
    self.value = [float(u) for u in value.rstrip().split()]

  def __repr__(self):
    """ Dumps a representation of self. """
    args = []
    if 'keyword' not in self.__class__.__dict__ and 'keyword' in self.__dict__:
      args.append('keyword={0.keyword!r}'.format(self))
    if 'shape' not in self.__class__.__dict__ and 'shape' in self.__dict__: 
      args.append('shape={0.shape!r}'.format(self))
    if 'units' not in self.__class__.__dict__ and 'units' in self.__dict__: 
      args.append('units={0.units!r}'.format(self))
    if len(getattr(self, 'shape', ())) > 0: 
      args.append('value={0!r}'.format(self.value.magnitude))
    else: args.append('value={0}'.format(float(self.value)))
    return '{0.__class__.__name__}({1})'.format(self, ', '.join(args))

