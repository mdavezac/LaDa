""" Defines decorators which can cache to file. """
__docformat__ = "restructuredtext en"

class FileCache(object):
  """ Caches results to a file in current directory.

      If the file exists, reloads and returns corresponding result.
      If the first argument contains a ``dirattr`` attribute ("directory" by
      default), then the file should in that directory.
  """
  def __init__(self, filename, dirattr="directory"):
    """ Initializes filecache objet. """
    self.filename = filename
    """ Filename where to cache results. """
    self.dirattr  = dirattr
    """ Name of the attribute with the directory where to cache file. """
  def __call__(self, method):
    """ Creates decorated function. """

    def wrapper(this, *args, **kwargs):
      """ Wrapper around cached method. """
      from os.path import exists
      from pickle import load, dump
      from ..mpi import Communicator

      # gesses at where the communicator is.
      comm = Communicator(kwargs.get('comm', getattr(this, 'comm', None)))
      # gets path to file, including directory.
      path = self._getpath(this)
      # if the file exists, then reloads result.
      if exists(path):
        calculation = None
        if comm.is_root: 
          with open(path, 'r') as file: calculation = load(file)
        calculation = comm.broadcast(calculation)
        return calculation

      # if the file does not exist, create result and cache it.
      calculation = method(this, *args, **kwargs)
      if calculation is None: return None
      if comm.is_root:
        with open(path, 'w') as file: dump(calculation, file)
      return calculation

    result = wrapper 
    result.__name__ = method.__name__
    result.__doc__ = method.__doc__
    result.__module__ = method.__module__
    result.__dict__.update(method.__dict__)
    return result

  def _getpath(self, this):
    """ Returns path to cache file. """
    from os.path import join
    if self.dirattr is not None:
      try: directory = getattr(this, self.dirattr)
      except: path = self.filename
      else:   path = join(directory, self.filename)
    return path

  def uncache(self, this):
    """ Removes file if it exists. """
    from os import remove
    from os.path import exists
    if exists(path): remove(path)

  
