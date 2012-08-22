from .keywords import BaseKeyword
class AttrBlock(BaseKeyword):
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
    self.__dict__['_input'] = {}

    # then call base constructor
    super(AttrBlock, self).__init__(keyword=keyword, raw=raw)

    # now so we get doctstrings right.
    self._input = {}
    """ Dictionary of crystal inputs. """
    

  def __getattr__(self, name):
    """ passes through the input keywords in :py:attr:`_input`. """
    from ..error import AttributeError
    if name not in self._input: 
      raise AttributeError('Unknown attribute {0}.'.format(name))
    result = self._input[name]
    return result.__get__(self) if hasattr(result, '__get__') else result
  def __setattr__(self, name, value):
    """ passes through the input keywords in :py:attr:`_input`. 
    
        If the input value is derived from
        :py:class:`~lada.tools.keyword.BaseKeyword`, then it is added to
        :py:attr:`_input`. Otherwise, super is called.
    """
    if isinstance(value, BaseKeyword): 
      self._input[name] = value
      if not hasattr(value, 'keyword'): self._input[name].keyword = name
    elif name in self._input:
      result = self._input[name]
      if hasattr(result, '__set__'): result.__set__(self, value)
      else: self._input[name] = value
    else: super(AttrBlock, self).__setattr__(name, value)
  def __delattr__(self, name):
    """ passes through the input keywords in :py:attr:`_input`.  """
    if name in self._input: del self._input[name]
    else: super(AttrBlock, self).__delattr__(name)
  def __dir__(self):
    """ List of attributes and members. """
    return list( set(self._input.iterkeys())                                   \
                 | set(self.__dict__.iterkeys())                               \
                 | set(dir(self.__class__)) )

  def add_keyword(self, name, value=None):
    """ Adds/Sets input keyword. """
    # if known keyword, then go through setattr mechanism.
    # this makes sure we recognize the type of value and the already registered
    # keyword.
    if name in self._input:  setattr(self, name, value)
    # if value is None, then transform it to True. 
    # This is a keyword which is either there or not there, like EXTPRT.
    elif value is None: self._input[name] = True
    # boolean case
    elif value is True or value is False:
      self._input[name] = value
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
        # if made up of string, then go back to string.
        if all(isinstance(u, str) for u in n): n = [lines]
        # if only one element use that rather than list
        if len(n) == 1: n = n[0]
      # if multiple line, keep as such
      else: n = value
      self._input[name] = n
    # otherwise, just set the keyword.
    else: self._input[name] = value
    # return self to allow chaining calls.
    return self

  def __repr__(self, defaults=False, name=None):
    """ Representation of this instance. """
    from .uirepr import uirepr
    defaults = self.__class__() if defaults else None
    return uirepr(self, name=name, defaults=defaults)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Creates user friendly representation. """
    from .uirepr import template_ui_repr, add_to_imports

    results = template_ui_repr(self, imports, name, defaults, exclude)
    if name is None:
      name = getattr(self, '__ui_name__', self.__class__.__name__.lower())

    for key, value in self._input.iteritems():
      if exclude is not None and key in exclude: continue
      if hasattr(value, '__ui_repr__'): 
        default = None if defaults is None                                     \
                  else defaults._input.get(key, None)
        newname = name + '.' + key
        partial = value.__ui_repr__(imports, newname, default)
        results.update(partial)
        if newname in results:              doinit = False
        elif default is None:               doinit = True
        else: doinit = type(value) is not type(default)
        if doinit:
          results[newname] = '{0.__class__.__name__}()'.format(value)
          add_to_imports(value, imports)
      elif isinstance(value, BaseKeyword):
        value = getattr(self, key)
        string = repr(value)
        if defaults is not None and key in defaults._input                     \
           and type(value) is type(getattr(defaults, key))                     \
           and string == repr(getattr(defaults, key)): continue
        key = '{0}.{1}'.format(name, key) 
        results[key] = string
        add_to_imports(value, imports)
      elif value is None:
        if defaults is not None and key in defaults._input                     \
           and defaults._input[key] is None: continue
        results['{0}.add_keyword({1!r})'.format(name, key)] = None
      else:
        if defaults is not None and key in defaults._input                     \
           and type(value) is type(defaults._input[key])                       \
           and repr(value) == repr(defaults._input[key]): continue
        results['{0}.add_keyword({1!r}, {2!r})'.format(name, key, value)]      \
            = None
        add_to_imports(value, imports)
    
    return results
  
  def __getstate__(self):
    d = self.__dict__.copy()
    crysinput = d.pop('_input')
    return d, crysinput
  def __setstate__(self, value):
    self.__dict__['_input'] = value[1]
    self.__dict__.update(value[0])

  def output_map(self, **kwargs):
    """ Map of keyword, value """
    result = {}
    for key, value in self._input.iteritems():
      if value is None: continue
      elif isinstance(value, bool):
        result[key] = '.TRUE.' if value else '.FALSE.'
      elif hasattr(value, 'output_map'):
        dummy = value.output_map(**kwargs)
        if dummy is not None: result.update(dummy)
      elif getattr(value, 'raw', None) is not None:
        result[key] = str(value.raw)
      elif hasattr(value, '__iter__'):
        result[key] =  ' '.join(str(u) for u in value)
      else: result[key] = str(value)
    return result

