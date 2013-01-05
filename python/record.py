""" Object with on-disk attributes. 

    Declares a class with attributes which are pickled on disk. 
"""
class Record(object):
  """ Object whose attributes are stored in a file. """
  def __init__(self, path=None):
    """ Initializes an on-disk object. """
    from .opt import RelativeDirectory
    if path is None: path = ".pylada_record"
    self.__dict__['path'] = RelativeDirectory(path).path
    """ Path to on-disk file. """

  def __getattr__(self, name):
    """ Reads on-disk attributes. """
    from pickle import load
    from os.path import exists
    from .opt import LockFile
    if not exists(self.path): raise AttributeError("Attribute {0} does not exist.".format(name))
    with LockFile(self.path) as lock:
      with open(self.path, 'r') as file:
        keys = load(file)
        if name not in keys: raise AttributeError("Attribute {0} does not exist.".format(name))
        return load(file)[1][name]

  def __setattr__(self, name, value):
    """ Writes on-disk attributes. """
    from pickle import dump, dumps, load
    from os.path import exists
    from .opt import LockFile
    try: dumps(value)
    except: raise RuntimeError("{0} is not a pickle-able object.".format(name))
    with LockFile(self.path) as lock:
      if exists(self.path): 
        with open(self.path, 'r') as file:
          dummy = load(file)
          dummy, dictionary = load(file)
      else: dictionary = {}
      dictionary[name] = value
      with open(self.path, 'w') as file:
        dump( set(dictionary.keys()), file )
        dump( ("This is a record.", dictionary), file)

  def __delattr__(self, name):
    """ Removes on-disk attribute. """
    from pickle import dump, load
    from os.path import exists
    from .opt import LockFile
    if not exists(self.path): return
    with LockFile(self.path) as lock:
      with open(self.path, 'r') as file:
        dummy = load(file)
        dummy, dictionary = load(file)
      if name not in dictionary: return 
      dictionary.pop(name)
      with open(self.path, 'w') as file:
        dump( set(dictionary.keys()), file )
        dump( ("This is a record.", dictionary), file)


  def __dir__(self):
    """ Returns list of on-disk attributes. """
    return [u for u in self._allattr() if u[0] != '_']

  def _allattr(self):
    """ Returns all attributes, both private and public. """
    from os.path import exists
    from pickle import load
    from .opt import LockFile
    with LockFile(self.path) as lock:
      if not exists(self.path): return []
      with open(self.path, 'r') as file: return list(load(file))
