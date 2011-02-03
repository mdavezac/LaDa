""" Mixin classes for extraction objects. """
__docformat__  = 'restructuredtext en'


class IOMixin(object):
  """ A mixin base clase which controls file IO. 

      Defines special property with file-like behaviors. 
      Makes it easier to change the behavior of the extraction class.
  """
  def __init__(self, directory=None, OUTCAR=None, FUNCCAR=None, CONTCAR=None):
    """ Initializes the extraction class. 

        :Parameters: 
          directory : str or None
            path to the directory where the VASP output is located. If none,
            will use current working directory. Can also be the path to the
            OUTCAR file itself. 
          OUTCAR : str or None
            If given, this name will be used, rather than files.OUTCAR.
          CONTCAR : str or None
            If given, this name will be used, rather than files.CONTCAR.
          FUNCCAR : str or None
            If given, this name will be used, rather than files.FUNCCAR.
    """
    from .. import files
    
    object.__init__(self)

    self.OUTCAR  = OUTCAR if OUTCAR != None else files.OUTCAR
    """ Filename of the OUTCAR file from VASP. """
    self.CONTCAR  = CONTCAR if CONTCAR != None else files.CONTCAR
    """ Filename of the CONTCAR file from VASP. """
    self.FUNCCAR  = FUNCCAR if FUNCCAR != None else files.FUNCCAR
    """ Filename of the FUNCCAR file containing the pickled functional. """

  def __outcar__(self):
    """ Returns path to OUTCAR file.

        :raise IOError: if the OUTCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.OUTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  def __funccar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.FUNCCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')

  def __contcar__(self):
    """ Returns path to FUNCCAR file.

        :raise IOError: if the FUNCCAR file does not exist. 
    """
    from os.path import exists, join
    path = join(self.directory, self.CONTCAR)
    if not exists(path): raise IOError("Path {0} does not exist.\n".format(path))
    return open(path, 'r')


class SearchMixin(object):
  """ Defines search methods. """
  def __init__(self): object.__init__(self)
  
  def _search_OUTCAR(self, regex):
    """ Looks for all matches. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    result = []
    regex  = compile(regex)
    with self.__outcar__() as file:
      for line in file: 
        found = regex.search(line)
        if found != None: yield found

  def _find_first_OUTCAR(self, regex):
    """ Returns first result from a regex. """
    for first in self._search_OUTCAR(regex): return first
    return None

  def _rsearch_OUTCAR(self, regex):
    """ Looks for all matches starting from the end. """
    from os.path import exists, join
    from re import compile
    from numpy import array

    result = []
    regex  = compile(regex)
    with self.__outcar__() as file: lines = file.readlines()
    for line in lines[::-1]:
      found = regex.search(line)
      if found != None: yield found

  def _find_last_OUTCAR(self, regex):
    """ Returns first result from a regex. """
    for last in self._rsearch_OUTCAR(regex): return last
    return None
