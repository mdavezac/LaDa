from .keywords import BaseKeyword

class ListBlock(BaseKeyword, list):
  """ Defines block input acting as a list.
  
      This type of group input can contain subitems arranged in a list (rather
      than accessible as attributes, as in
      :py:class:`lada.tools.input.block.AttrBlock`)
  """
  __ui_name__ = 'structure'
  """ Name used in user-friendly representation. """
  def __init__(self, keyword=None, raw=None):
    """ Creates a block. 

        :param str keyword:
          Keyword indicating the name of the block.
    """
    BaseKeyword.__init__(self, keyword=keyword, raw=raw)
    list.__init__(self)

  def append(self, keyword, raw=None):
    """ Appends an item.

        :return: self, so calls can be chained. 
    """
    from ...error import ValueError
    if isinstance(keyword, str):
      keyword = BaseKeyword(keyword=keyword, raw=raw)
    elif not isinstance(keyword, BaseKeyword):
      raise ValueError('Wrong argument to {0.keyword}.append().'.format(self))
    list.append(self, keyword)
    return self

  
  def output_map(self, **kwargs):
    """ Print input to crystal. """
    from .tree import Tree
    root = Tree()
    result = root if getattr(self, 'keyword', None) is None                    \
             else root.descend(self.keyword)
    if getattr(self, 'raw', None) is not None:
      result.prefix = str(self.raw)
    for item in self:
      keyword = getattr(item, 'keyword', None)
      value   = getattr(item, 'raw', None)
      if hasattr(item, 'output_map'):
        dummy = item.output_map(**kwargs)
        if dummy is not None: result.update(dummy)
      elif value is not None:
        result[keyword] = value
      elif hasattr(value, '__iter__'):
        result[keyword] =  ' '.join(str(u) for u in value)
      else: result[keyword] = str(value)
    return root

  def read_input(self, tree, owner=None, **kwargs):
    """ Sets object from input tree. """
    from ...error import IndexError
    from .tree import Tree
    for key, value in tree.iteritems():
      key = key.lower()
      try: 
        if isinstance(value, Tree) and key not in self:
            raise IndexError('Unknown group {0}'.format(key))
        self.append(key, value)
      except:
        from sys import exc_info
        type, value, traceback = exc_info()
        message = 'ERROR when reading {0}.'.format(key)
        if value is None: type.args = type.args, message
        else: value = value, message
        raise type, value, traceback
  

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ Dumps representation to string. """
    from ..uirepr import add_to_imports

    result = super(ListBlock, self).__repr__()
    indent = ' '.join('' for i in xrange(result.find('(')+1))
    add_to_imports(self, imports)

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

  def __repr__(self, indent=''): 
    """ Dumps representation to string. """
    return self.__ui_repr__({}).itervalues().next()
