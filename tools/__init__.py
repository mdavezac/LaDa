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
      for key in kwargs.keys():
        if key in ignore: continue
        # direct attributes.
        if hasattr(self, key): setattr(self, key, kwargs.pop(key))
        else:
          found = False
          for other in setothers:
            if hasattr(self, other) and hasattr(getattr(self, other), key):
              setattr(getattr(self, other), key, kwargs.pop(key))
              found = True
              break
          if found == False:
            raise ValueError( "Unkwown keyword argument to {0.__class__.__name__}: {1}"\
                              .format(self, key) )
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


class SuperCall(object):
  """ Obviates issues when using a "super" functional.

      Since functionals of a job-folder are deepcopied, the following line
      will not result in calling the next class in the __mro__.
  
      >>> jobfolder.functional = super(Functional, functional)

      Indeed, this line will first call the __getitem__, __setitem__ (or
      __deepcopy__) of the super object. In general, this means we end-up with
      ``jobfolder.function == functional``.

      This class obviates this difficulty.

      >>> jobfolder.functional = SuperCall(Functional, functional)
  """
  def __init__(self, class_, object_):
    object.__init__(self)
    self.__dict__['_class'] = class_
    self.__dict__['_object'] = object_
  def __call__(self, *args, **kwargs):
    return super(self._class, self._object).__call__(*args, **kwargs)
  def __getattr__(self, name):
    try: return getattr(super(self._class, self._object), name)
    except: return getattr(self._object, name)
  def __setattr__(self, name, value): setattr(self._object, name, value)
  def __getstate__(self): return self._class, self._object
  def __setstate__(self, args):
    self.__dict__['_class'] = args[0]
    self.__dict__['_object'] = args[1]
  def __dir__(self): return dir(super(self._class, self._object))
  def __repr__(self): return repr(super(self._class, self._object))
  def copy(self):
    """ Performs deepcopy of self. """
    from copy import deepcopy
    class_ = deepcopy(self.__dict__['_class'])
    object_= deepcopy(self.__dict__['_object'])
    return self.__class__(class_, object_)
