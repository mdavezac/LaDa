""" Holds base classes and mixins for extraction objects. """
__docformat__ = "restructuredtext en"
__all__ = ['AbstractExtractBase', 'search_factory']
from abc import ABCMeta, abstractproperty

class AbstractExtractBase(object):
  """ Abstract base class for extraction classes. 
  
      Defines a number of members common to all extraction classes:
        - directory: root directory where output should exist.
  """
  __metaclass__ = ABCMeta
  def __init__(self, directory=None):
    """ Initializes an extraction base class.

        :Parameters: 
          directory : str or None
            Root directory for extraction. If None, will use current working directory.
    """
    object.__init__(self)

    from os import getcwd
    from . import RelativeDirectory

    if directory is None: directory = getcwd()
    self._directory = RelativeDirectory(directory, hook=self.__directory_hook__)
    """ Directory where output should be found. """

  @property
  def directory(self):
    """ Directory where output should be found. """
    return self._directory.path
  @directory.setter
  def directory(self, value): self._directory.path = value

  @abstractproperty
  def success(self):
    """ Checks for success. 

        Should never ever throw!
        True if calculations were successfull, false otherwise.
    """
    pass


  def __directory_hook__(self):
    """ Called whenever the directory changes. """
    self.uncache()

  def uncache(self): 
    """ Uncache values. """
    self.__dict__.pop("_properties_cache", None)

  def __copy__(self):
    """ Returns a shallow copy of this object. """
    from . import RelativeDirectory
    result = self.__class__(directory=self.directory)
    result.__dict__ = self.__dict__.copy()
    result._directory = RelativeDirectory( self._directory.path,\
                                           self._directory._envvar, 
                                           result.uncache )
    return result

  def copy(self, **kwargs):
    """ Returns a shallow copy of this object.

        :param kwargs:
          Any keyword argument is set as an attribute of this object.
          The attribute must exist.
    """
    result = self.__copy__()
    for k, v in kwargs.iteritems():
      if not hasattr(self, k): raise RuntimeError('Attribute {0} does not exist.'.format(k))
      setattr(result, k, v)
    return result

  def __getstate__(self):
    d = self.__dict__.copy()
    if "_directory" in d: d["_directory"].hook = None
    return d

  def __setstate__(self, arg):
    self.__dict__.update(arg)
    if hasattr(self, "_directory"): self._directory.hook = self.uncache

  def __repr__(self):
    return "{0}(\"{1}\")".format(self.__class__.__name__, self._directory.unexpanded)

def search_factory(name, filename, module):
  """ Factory to create Mixing classes capable of search a given file. """
  doc = \
    """ A mixin to include standard methods to search {0}.
    
        This mixin only includes the methods themselves. It expects the derived
        class to have an {0} attribute. 
    """.format(filename.upper())
  def __outcar__(self):
    """ Returns path to OUTCAR file.

        :raise IOError: if the OUTCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, getattr(self, filename.upper()))
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')
  __outcar__.__name__ = '__{0}__'.format(filename.lower())

  def _search_OUTCAR(self, regex, flags=0):
    """ Looks for all matches. """
    from re import compile, M as moultline

    regex  = compile(regex, flags)
    with getattr(self, __outcar__.__name__)() as file:
      if moultline & flags: 
        for found in regex.finditer(file.read()): yield found
      else:
        for line in file: 
          found = regex.search(line)
          if found is not None: yield found
  _search_OUTCAR.__name__ = '_search_{0}'.format(filename.upper())

  def _find_first_OUTCAR(self, regex, flags=0):
    """ Returns first result from a regex. """
    for first in getattr(self, _search_OUTCAR.__name__)(regex, flags): return first
    return None
  _find_first_OUTCAR.__name__ = '_find_first_{0}'.format(filename.upper())

  def _rsearch_OUTCAR(self, regex, flags=0):
    """ Looks for all matches starting from the end. """
    from re import compile, M as moultline

    regex  = compile(regex)
    with getattr(self, __outcar__.__name__)() as file:
      lines = file.read() if moultline & flags else file.readlines()
    if moultline & flags: 
      for v in [u for u in regex.finditer(lines)][::-1]: yield v
    else:
      for line in lines[::-1]:
        found = regex.search(line)
        if found is not None: yield found
  _rsearch_OUTCAR.__name__ = '_rsearch_{0}'.format(filename.upper())

  def _find_last_OUTCAR(self, regex, flags=0):
    """ Returns first result from a regex. """
    for last in getattr(self, _rsearch_OUTCAR.__name__)(regex, flags): return last
    return None

  attrs = { __outcar__.__name__: __outcar__,
            _search_OUTCAR.__name__: _search_OUTCAR,
            _rsearch_OUTCAR.__name__: _rsearch_OUTCAR,
            _find_first_OUTCAR.__name__: _find_first_OUTCAR,
            _find_last_OUTCAR.__name__: _find_last_OUTCAR,\
            '__doc__': doc,
            '__module__': module }
  return type(name, (), attrs)
