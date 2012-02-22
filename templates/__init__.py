""" Miscellaneous ressources to create functionals. """
from collections import namedtuple


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
    from ..opt import RelativeDirectory
    structure = deepcopy(structure)
    self      = deepcopy(self)
    outdir    = getcwd() if outdir is None else RelativeDirectory(outdir).path
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
      extract = self.Extract(comm = comm, outcar = outdir)
      if extract.success: return extract # in which case, returns extraction object.
    return function(self, *wargs, **kwargs)
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
  if hasattr(self, '_properties_cache'): del ob._properties_cache


Program = namedtuple('Program', ['program', 'cmdline', 'directory', 'stdout', 'stderr'])
""" Holds data about program to launch. """

def execute_program(program, append=False, **kwargs):
   """ Executes a program. """
   from subprocess import Popen
   from shlex import split as split_cmd
   from lada import mpirun_exe as mpirun_exe_global
   from lada.opt import ChangeDirectory
   
   d = kwargs.copy()
   d['program'] = program
   d['cmdline'] = cmdline
   if mpirun_exe == None: mpirun_exe = mpirun_exe_global
   cmd = mpirun_exe.format(d)
   
   with ChangeDirectory(program.directory) as cwd:
     file_out = open(program.stdout, "a" if append else "w") if out is not None else None 
     file_err = open(program.stderr, "a" if append else "w") if err is not None else None 
     try:
       vasp_proc = Popen(split_cmd(cmd), stdout=file_out, stderr=file_err)
       vasp_proc.wait()
     finally:
       file_out.close()
       file_err.close()

def add_setter(method, docstring = None): 
  """ Adds an input-like setter property. """
  def _not_available(self): raise RuntimeError("Error: No cheese available.")
  if docstring is None and hasattr(method, "__doc__"): docstring = method.__doc__
  return property(fget=_not_available, fset=method,  doc=docstring)
