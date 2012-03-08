""" Miscellaneous ressources to create functionals. """
def stateless(function):
  """ Decorator to make a function stateless.
  
      Deepcopies structure and self parameters.
      Also sets outdir to getcwd if it is None.
  """
  from functools import wraps

  @wraps(function)
  def wrapper(self, structure, outdir=None, **kwargs ):
    from copy import deepcopy
    from os import getcwd
    from ..misc import RelativePath
    structure = deepcopy(structure)
    self      = deepcopy(self)
    outdir    = getcwd() if outdir is None else RelativePath(outdir).path
    return function(self, structure, outdir, **kwargs)
  return wrapper

def assign_attributes(setothers=None, ignore=None):
  """ Decorator to assign keywords to attributes. """
  from functools import wraps
  if setothers == None: setothers = []
  if ignore == None: ignore = set()

  @wraps(assign_attributes)
  def decorator(function):
    """ Decorator to assign keywords to attributes. """
    @wraps(function)
    def wrapper(self, structure, outdir=None, comm=None, **kwargs):
      # if other keyword arguments are present, then they are assumed to be
      # attributes of self, with value to be changed before launch. 
      for key, value in kwargs.items():
        if key in ignore: continue
        # direct attributes.
        if hasattr(self, key): setattr(self, key, value)
        # properties attributes.
        elif hasattr(self.__class__, key): setattr(self, key, value)
        else:
          found = False
          for other in setothers:
            if hasattr(self, other) and hasattr(getattr(self, other), key):
              setattr(getattr(self, other), key, value)
              found = True
              break
          if found == False:
            raise ValueError( "Unkwown keyword argument to {0.__class__.__name__}: {1}={2}"\
                              .format(self, key, value) )
      return function(self, structure, outdir=outdir, comm=comm, **kwargs)
    return wrapper
  return decorator

def check_success(function):
  """ Decorator to check for success prior to running functional. """
  from functools import wraps
  @wraps(function)
  def wrapper(self, *args, **kwargs):
    # Checks for previous run, or deletes previous run if requested.
    if not kwargs.get('overwrite', False):
      extract = self.Extract(outcar = kwargs['outdir'])
      if extract.success: return extract # in which case, returns extraction object.
    return function(self, *args, **kwargs)
  return wrapper

def check_success_generator(function):
  """ Decorator to check for success prior to running functional. 
  
      Generator version. Yields stuff.
  """
  from functools import wraps
  @wraps(function)
  def wrapper(self, *args, **kwargs):
    # Checks for previous run, or deletes previous run if requested.
    if not kwargs.get('overwrite', False):
      extract = self.Extract(outcar = kwargs['outdir'])
      if extract.success: yield extract # in which case, returns extraction object.
      return 
    for n in function(self, *args, **kwargs): yield n
  return wrapper

def make_cached(method):
  """ Caches the result of a method for futur calls. """
  from functools import wraps

  @wraps(method)
  def wrapped(*args, **kwargs):
    if not hasattr(args[0], '_properties_cache'): 
      setattr(args[0], '_properties_cache', {}) 
    cache = getattr(args[0], '_properties_cache')
    if method.__name__ not in cache:
      cache[method.__name__] = method(*args, **kwargs)
    return cache[method.__name__]
  return wrapped

def uncache(ob):
  """ Uncaches results cached by @make_cached. """ 
  if hasattr(ob, '_properties_cache'): del ob._properties_cache


