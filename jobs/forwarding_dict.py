""" Special dictioanry class to forward attributes indefinitely. 

    In practice, this class looks like a namespace which returns further namespaces. 
"""
__docformat__ = "restructuredtext en"
__all__ = ['ForwardingDict']
from collections import MutableMapping

class ForwardingDict(MutableMapping): 
  """ An *ordered* dictionary which forwards attributes and calls. """
  def __init__(self, dictionary = None, _attr_list=None, ordered=True, **kwargs):
    """ Initializes a ForwardingDict instance. """
    from ..opt import OrderedDict
    from . import default_params as default
    self._is_initializing_forwarding_dict = True
    """ Tells get/setattr that Forwarding dict is being initialized. """
    super(ForwardingDict, self).__init__()

    self.readonly      = kwargs.pop('readonly', default.readonly)
    """ Whether items can be modified in parallel using attribute syntax. """
    self.naked_end     = kwargs.pop('naked_end', default.naked_end)
    """ Whether last item is returned as is or wrapped in ForwardingDict. """
    self.only_existing = kwargs.pop('only_existing', default.only_existing)
    """ Whether attributes can be added or only modified. """
    self._attr_list    = [] if _attr_list == None else _attr_list
    """ List of attributes of attributes, from oldest parent to youngest grandkid. """
    dicttype = OrderedDict if ordered else dict
    self.dictionary    = dicttype({} if dictionary == None else dictionary)
    """" The dictionary for which to unroll attributes. """
    del self._is_initializing_forwarding_dict
    assert len(kwargs) == 0, ValueError("Unkwnown keyword arguments:{0}.".format(kwargs.keys()))

  def call(self, *args, **kwargs):
    """ Forwards call to items in dictionary.
    
        In practice, only forward calls to those items which can be completely
        unrolled and for which the unrolled object is callable.

        There seems to be a bug in python 2.6 or 2.7 which make instances
        derived from object with a __call__ bound method in __dict__ (rather
        than in __class__.__dict__) uncallable. Linked to tp_getattr deprecation?
    """
    result = self.copy(attr_list=[])
    result.clear()
    for key, value in self.dictionary.iteritems():
      if hasattr(value, "__call__"): result[key] = value(*args, **kwargs)
    if self.naked_end and len(result) == 1: return result[result.keys()[0]]
    return result


  @property
  def parent(self):
    """ Returns a ForwardingDict with parent items of self, eg unrolled once. """
    return self.copy(_attr_list=self._attr_list[:-1])

  @property
  def root(self):
    """ Returns a ForwardingDict with root grandparent. """
    return self.copy(_attr_list=[])

  @property
  def _attributes(self):
    """ Returns attributes special to this ForwardingDict. """
    from functools import reduce
    from itertools import chain

    result = set()
    attrs = len(self._attr_list) > 0
    for value in self.dictionary.itervalues():
      if attrs: value = reduce(getattr, chain([value], self._attr_list))
      result |= set(dir(value))
      
    return result

  def __getattr__(self, name):
    """ Returns a Forwarding dict with next requested attribute. """
    from functools import reduce
    from itertools import chain

    if name not in self._attributes:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )
    attrs = len(self._attr_list) > 0
    result = self.copy(append=name)
    for key, value in self.dictionary.iteritems():
      if attrs: value = reduce(getattr, chain([value], self._attr_list))
      if not hasattr(value, name): del result[key]
    if self.naked_end and len(result.dictionary) == 1: return result[result.keys()[0]]
    if len(result.dictionary) == 0: 
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )
    return result

  def __setattr__(self, name, value):
    """ Forwards attribute setting. """
    from functools import reduce
    from itertools import chain
    # prior to initialization.
    if name == "_is_initializing_forwarding_dict": 
      super(ForwardingDict, self).__setattr__(name, value)
    if "_is_initializing_forwarding_dict" in self.__dict__:
      super(ForwardingDict, self).__setattr__(name, value)

    # After initialization.
    # First checks for attribute in ForwardingDict instance.
    try: super(ForwardingDict, self).__getattribute__(name)
    except AttributeError: pass
    else: super(ForwardingDict, self).__setattr__(name, value); return

    # checks this dictionary is writable.
    if self.readonly: raise RuntimeError("ForwardingDict instance is read-only.")

    # Case with no attributes to unroll.
    found = False
    attrs = len(self._attr_list) > 0
    for item in self.dictionary.values():
      if attrs: # unroll attribute list.
        try: item = reduce(getattr, chain([item], self._attr_list))
        except AttributeError: continue
      if hasattr(item, name) or not self.only_existing:
        found = True
        setattr(item, name, value)
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  def __delattr__(self, name):
    """ Deletes an attribute or forwarded attribute. """
    from functools import reduce
    from itertools import chain

    try: super(ForwardingDict, self).__delattr__(name)
    except AttributeError: pass
    else: return
    
    if self.readonly: raise RuntimeError("ForwardingDict instance is read-only.")

    found = False
    attrs = len(self._attr_list) > 0
    for item in self.dictionary.values():
      if attrs:
        try: item = reduce(getattr, chain([item], self._attr_list))
        except AttributeError: continue
      if hasattr(item, name): 
        delattr(item, name)
        found = True
    if not found:
      raise AttributeError( "Attribute {0} not found in {1} instance."\
                            .format(name, self.__class__.__name__) )

  def __dir__(self):
    from itertools import chain
    results = chain( [u for u in self.__dict__ if u[0] != '_'], \
                     [u for u in dir(self.__class__) if u[0] != '_'], \
                     self._attributes )
    return list(set(results))


  def __getitem__(self, key):
    from functools import reduce
    from itertools import chain
    if len(self._attr_list) == 0: return self.dictionary[key]
    return reduce(getattr, chain([self.dictionary[key]], self._attr_list))
  def __setitem__(self, key, value):
    """ Add/modify item to dictionary.

        Items can be truly added only to root dictionary.
    """
    from functools import reduce
    from itertools import chain
    # root dictioanary.
    if len(self._attr_list) == 0: self.dictionary[key] = value; return
    # checks this is writable.
    assert not self.readonly, RuntimeError("This ForwardingDict is readonly.")
    assert key in self.dictionary,\
           KeyError( "{0} is not in the ForwaringDict. Items "\
                      "cannot be added to a non-root ForwardingDict.".format(key))
    # non-root dict: must set innermost attribute.
    o = self.dictionary[key]
    if len(self._attr_list) > 1: 
      try: o = reduce(getattr, chain([o], self._attr_list[:-1]))
      except AttributeError:
        raise AttributeError( "Could not unroll list of attributes for object in {0}: {1}."\
                              .format(key, self._attr_list) )  
    assert (not self.only_existing) or hasattr(o, self._attr_list[-1]), \
           KeyError( "{0} cannot be set with current attribute list.\n{1}\n"\
                      .format(key, self._attr_list) )
    setattr(o, self._attr_list[-1], value)
  def __delitem__(self, key): 
    """ Removes item from dictionary. """
    o = self.dictionary[key]
    del self.dictionary[key]
    return o
  def __len__(self): return len(self.dictionary)
  def __contains__(self, key): return key in self.dictionary
  def __iter__(self): return self.dictionary.__iter__()
  def keys(self): return self.dictionary.keys()


  def __copy__(self):
    """ Returns a shallow copy of this object. """
    result = self.__class__()
    result.__dict__.update(self.__dict__)
    result.dictionary = self.dictionary.copy()
    return result

  def copy(self, append=None, dict=None, **kwargs):
    """ Returns a shallow copy of this object.
     
        :Parameters:
          append : str or None
            If not none, will append value to list attributes of the copy. In
            that case, the list of attributes ``_attr_list`` is copied
            *deeply*.
          kwargs : dict
            Any other attribute to set in the ForwardingDict instance. Note
            that only attributes of the ForwardingDict instance are
            set/modified. This is npt propagated to the object the dict holds.
    """
    from copy import copy, deepcopy
    result = copy(self)
    assert append == None or "_attr_list" not in kwargs,\
           ValueError( "Cannot copy attribute _attr_list as "\
                        "a keyword and as ``append`` simultaneously." )
    if 'dictionary' in kwargs: result.dictionary = kwargs.pop('dictionary').copy()
    for key, value in kwargs.iteritems():
      super(ForwardingDict, result).__setattr__(key, value)

    if append != None:
      result._attr_list = deepcopy(self._attr_list)
      result._attr_list.append(append)

    return result


  def __str__(self): return self.dictionary.__str__()
  def __repr__(self): 
    """ Returns a represetation of this object. """
    if len(self) == 0: return '{}'
    if len(self) == 1: return "{{'{0}': {1}}}".format(self.keys()[0], repr(self.values()[0]))
    string = "{\n"
    m = max(len(k) for k in self.keys())
    for k, v in self.iteritems():
      string += "  '{0}': {2}{1},\n".format(k, repr(v), "".join(" " for i in range(m-len(k))))
    return string + "}"

  def __setstate__(self, state):
    """ Reloads the state from a pickle. 

        This is defined explicitely since otherwise, the call would go through
        __getattr__ and start an infinite loop.
    """ 
    self.__dict__.update(state)
