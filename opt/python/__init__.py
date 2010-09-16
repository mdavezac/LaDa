""" Miscellaneous """
import _opt 
from contextlib import contextmanager
__load_vasp_in_global_namespace__ = _opt.__load_vasp_in_global_namespace__
__load_escan_in_global_namespace__ = _opt.__load_escan_in_global_namespace__
cReals = _opt.cReals
_RedirectFortran = _opt._RedirectFortran
streams = _opt._RedirectFortran.fortran
ConvexHull = _opt.ConvexHull
ErrorTuple = _opt.ErrorTuple

class _RedirectPy:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    self.append = append
  def __enter__(self):
    import sys
    if self.unit == streams.input:    self.old = sys.stdin
    elif self.unit == streams.error:  self.old = sys.stderr
    elif self.unit == streams.output: self.old = sys.stdout
    else: raise RuntimeError("Unknown redirection unit.")
    self.file = open(self.filename if len(self.filename)\
                else "/dev/null", "a" if self.append else "w")
    if self.unit == streams.input:    sys.stdin  = self.file
    elif self.unit == streams.error:  sys.stderr = self.file
    elif self.unit == streams.output: sys.stdout = self.file
    else: raise RuntimeError("Unknown redirection unit.")
    return self
  def __exit__(self, *wargs):
    import sys 
    if self.unit == streams.input:    sys.stdin  = self.old 
    elif self.unit == streams.error:  sys.stderr = self.old 
    elif self.unit == streams.output: sys.stdout = self.old 
    else: raise RuntimeError("Unknown redirection unit.")
    self.file.close()
      
class _RedirectAll:
  """ Redirects C/C++ input, output, error. """
  def __init__(self, unit, filename, append):
    self.unit = unit
    self.filename = filename
    if self.filename == None: self.filename = ""
    if self.filename == "/dev/null": self.filename = ""
    self.append = append
  def __enter__(self):
    import os
    import sys
    if self.unit == streams.input:    self.old = os.dup(sys.__stdin__.fileno())
    elif self.unit == streams.output: self.old = os.dup(sys.__stdout__.fileno())
    elif self.unit == streams.error:  self.old = os.dup(sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")
    if len(self.filename) == 0:
      self.filefd = os.open(os.devnull, os.O_RDWR | os.O_APPEND)
    else:
      self.filefd = os.open(self.filename, (os.O_RDWR | os.O_APPEND) if self.append\
                                           else (os.O_CREAT | os.O_WRONLY | os.O_EXCL))
    if self.append: self.file = os.fdopen(self.filefd, 'a')
    else: self.file = os.fdopen(self.filefd, 'w')
    if self.unit == streams.input:    os.dup2(self.file.fileno(), sys.__stdin__.fileno())
    elif self.unit == streams.output: os.dup2(self.file.fileno(), sys.__stdout__.fileno())
    elif self.unit == streams.error:  os.dup2(self.file.fileno(), sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")
    return self
  def __exit__(self, *wargs):
    import os
    import sys
    try: self.file.close()
    except: pass
    if self.unit == streams.input:    os.dup2(self.old, sys.__stdin__.fileno()) 
    elif self.unit == streams.output: os.dup2(self.old, sys.__stdout__.fileno())
    elif self.unit == streams.error:  os.dup2(self.old, sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")

def redirect_all(output=None, error=None, input=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      :Parameters:
        output
          Filename to which to redirect output. 
        error
          Filename to which to redirect err. 
        input
          Filename to which to redirect input. 
        append
          If true, will append to files. All or nothing.
  """
  from contextlib import nested
  result = []
  for value, unit in [ (output, streams.output), (error, streams.error), (input, streams.input) ]:
    if value == None: continue
    result.append( _RedirectAll(unit=unit, filename=value, append=append) )
  return nested(*result)


def redirect(fout=None, ferr=None, fin=None, pyout=None, pyerr=None, pyin=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      :Parameters:
        fout
          Filename to which to redirect fortran output. 
        ferr
          Filename to which to redirect fortran err. 
        fin
          Filename to which to redirect fortran input. 
        pyout
          Filename to which to redirect C/C++ output. 
        pyerr
          Filename to which to redirect C/C++ err. 
        pyin
          Filename to which to redirect C/C++ input. 
        append
          If true, will append to files. All or nothing.
  """
  from contextlib import nested
  result = []
  for value, unit in [ (fout, streams.output), (ferr, streams.error), (fin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectFortran(unit=unit, filename=value, append=append) )
  for value, unit in [ (pyout, streams.output), (pyerr, streams.error), (pyin, streams.input) ]:
    if value == None: continue
    result.append( _RedirectPy(unit=unit, filename=value, append=append) )
  return nested(*result)



def read_input(filename, global_dict=None, local_dict = None, paths=None, comm = None):
  """ Executes input script and returns local dictionary (as class instance). """
  # stuff to import into script.
  from os import environ
  from os.path import abspath, expanduser, join
  from math import pi 
  from numpy import array, matrix, dot, sqrt, abs, ceil
  from numpy.linalg import norm, det
  from lada.crystal import Lattice, Site, Atom, Structure, fill_structure, FreezeCell, FreezeAtom
  from lada import physics
  from boost.mpi import world
  
  # Add some names to execution environment.
  if global_dict == None: global_dict = {}
  global_dict.update( { "environ": environ, "pi": pi, "array": array, "matrix": matrix, "dot": dot,\
                        "norm": norm, "sqrt": sqrt, "ceil": ceil, "abs": abs, "Lattice": Lattice, \
                        "Structure": Structure, "Atom": Atom, "Site": Site, "physics": physics,\
                        "fill_structure": fill_structure, "world": world, "FreezeCell": FreezeCell, \
                        "FreezeAtom": FreezeAtom, "join": join, "abspath": abspath, \
                        "expanduser": expanduser})
  if local_dict == None: local_dict = {}
  # Executes input script.
  execfile(filename, global_dict, local_dict)

  # Makes sure expected paths are absolute.
  if paths != None:
    for path in paths:
      if path not in local_dict: continue
      local_dict[path] = abspath(expanduser(local_dict[path]))
    
  # Fake class which will be updated with the local dictionary.
  class Input(physics.__class__): 
    def __getattr__(self, name):
      raise AttributeError( "All out of cheese!\n"
                            "Required input parameter \"{0}\" not found in {1}." \
                            .format(name, self.__name__) )
    def __delattr__(self, name): raise RuntimeError("Cannot delete object from input namespace.")
    def __setattr__(self, name, value): raise RuntimeError("Cannot set/change object in input namespace.")
  result = Input(filename)
  result.__dict__.update(local_dict)
  return result


  
class LockFile(object):
  """ Gets an advisory lock for file C{filename}.

      Creates a lock directory named after a file. Relies on the presumably
      atomic nature of creating a directory. 
      *Beware* of mpi problems! L{LockFile} is (purposefully) not mpi aware.
      If used unwisely, processes will lock each other out.
  """
  def __init__(self, filename, timeout = None, sleep = 0.5):
    """ Creates a lock object. 

        :Parameters:
          timeout
            will raise a RuntimeError when calling `lock` if
            the lock could not be aquired within this time.
          sleep
            Time to sleep between checks when waiting to acquire lock. 

        Does not acquire lock at this stage. 
    """
    from os.path import abspath, dirname, join
    self.filename = abspath(filename)
    """ Name of file to lock. """
    self.timeout = timeout
    """ Maximum amount of time to wait when acquiring lock. """
    self.sleep = sleep
    """ Sleep time between checks on lock. """
    self._owns_lock = False
    """ True if this object owns the lock. """

  def lock(self):
    """ Waits until lock is acquired. """
    from os import makedirs, error, mkdir
    from os.path import exists
    import time

    # creates parent directory first, if necessary.
    if not exists(self._parent_directory):
      try: makedirs(self._parent_directory) 
      except error: pass
    start_time = time.time()
    # loops until acqires lock.
    while self._owns_lock == False: 
      # tries to create director.
      try:
        self._owns_lock = True
        mkdir(self.lock_directory)
      # if fails, then another process already created it. Just keep looping.
      except error: 
        self._owns_lock = False
        if self.timeout != None:
          assert time.time() - start_time < self.timeout, \
                 RuntimeError("Could not acquire lock on %s." % (filename))
        time.sleep(self.sleep)

  def __enter__(self):
    """ Enters context. """
    assert self.timeout == None, RuntimeError("Cannot use LockFile as a context with timeout.")
    self.lock()
    return self

  def __exit__(self, *args):
    """ Exits context. """
    self.release()

  @property
  def lock_directory(self):
    """ Name of lock directory. """
    from os.path import join, basename
    return join(self._parent_directory, "." + basename(self.filename) + "-lada_lockdir")
 
  @property
  def _parent_directory(self):
    from os.path import abspath, dirname, join, basename
    return dirname(abspath(self.filename))

  @property
  def is_locked(self):
    """ True if a lock for this file exists. """
    from os.path import exists
    return exists(self.lock_directory)

  @property
  def owns_lock(self): 
    """ True if this object owns the lock. """
    return self._owns_lock

  def release(self):
    """ Releases a lock.

        It is an error to release a lock not owned by this object.
        It is also an error to release a lock which is not locked.
        Makes sure things are syncro. The latter is an internal bug though.
    """
    from os import rmdir
    assert self._owns_lock, IOError("Filelock object does not own lock.")
    assert self.is_locked, IOError("Filelock object owns an unlocked lock.")
    self._owns_lock = False
    rmdir(self.lock_directory)

  def __del__(self):
    """ Releases lock if still held. """
    if self.owns_lock and self.is_locked: self.release()

def acquire_lock(filename, sleep=0.5):
  """ Alias for a L{LockFile} context. 

      *Beware* of mpi problems! L{LockFile} is (purposefully) not mpi aware.
      Only the root node should use this method.
  """
  return LockFile(filename, sleep=sleep)

@contextmanager
def open_exclusive(filename, mode="r", sleep = 0.5):
  """ Opens file while checking for advisory lock.

      This context uses L{LockFile} to first obtain a lock.
      *Beware* of mpi problems! L{LockFile} is (purposefully) not mpi aware.
      Only the root node should use this method.
  """
  # Enter context.
  with LockFile(filename, sleep=sleep) as lock:
    yield open(filename, mode)

class RelativeDirectory(object):
  """ Directory property which is relative to the user's home.
  
      The directory wich is returned (eg __get__) is always absolute. However,
      it is stored relative to the user's home, and hence can be passed from
      one computer system to the next.
  """
  def __init__(self, name):
    """ Initializes the property. 
    
        :Param name: 
          Name of the property in the object.
    """
    import weakref

    self.im_owner = weakref.ref(None)
    self._value = "" if value == None else value
    """ Relative directory. """
    self.envvar = envvar if envvar != None else "~/"
    """ Fixed point. """
    self.hook = hook
    """ An object to call when the directory is changed.
    
        This callable need only take two argument, the 'self' of the owner and
        the absolute path to the directory.
    """

  @property.getter
  def envvar(self):
    """ Fixed point. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_envvar".format(self.name)
    if not hasattr(owner, name): return "~/"
    return getattr(owner, name)
  @property.setter
  def envvar(self, value):
    """ Fixed point. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_envvar".format(self.name)
    if value == None:
      if hasattr(owner, name): delattr(owner, name)
      return "~/"
    setattr(owner, name, value)

  @property.getter
  def hook(self):
    """ Call-back hook. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_hook".format(self.name)
    if not hasattr(owner, name): return lambda x: None
    return getattr(owner, name)
  @property.setter
  def hook(self, value):
    """ Call-back hook. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_hook".format(self.name)
    if value == None:
      if hasattr(owner, name): delattr(owner, name)
      return
    setattr(owner, name, value)

  @property.getter
  def directory(self):
    """ Directory relative to fixed point. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_directory".format(self.name)
    if not hasattr(owner, name): return lambda x: None
    return getattr(owner, name)
  @property.setter
  def directory(self, value):
    """ Directory relative to fixed point. """
    owner = self.im_owner()
    assert owner is None, RuntimeError("Owner is not known.")
    name = "_RelativeDirectory_{0}_directory".format(self.name)
    if value == None:
      if hasattr(owner, name):
        delattr(owner, name)
        self.hook(self.__get__(instance))
      return
    setattr(owner, name, value)
    self.hook(self.__get__(instance))
    

  def __repr__(self):
    """ Representation of relative directory. """
    owner = self.im_owner()
    if owner is None: return "{0}({\"1\"})".format(self.__class__.__name__, self.name)
    if self.envvar == "~/": return repr(self.directory)
    return "\"{0}\", \"{1}\"".format(self.directory, self.envvar)

  def __get__(self, instance, owner):
    """ Returns absolute directory """
    from os.path import expanduser, expandvars, join, normpath
    from weakref import ref
    
    if instance == None: return self
    self.im_owner = ref(instance)
    fxpoint = expandvars(expanduser(self.envvar(owner)))
    return normpath(join(fxpoint, self.directory))

  def __set__(self, instance, value):
    """ Sets directory relative to fixed-point. """
    from os.path import relpath
    from weakref import ref

    if instance == None: return
    self.im_owner ref(instance)
    # checking for a 3-tuple to set hook.
    if value != None and len(value) == 3 and hasattr(value[2], "__call__"):
      self.hook = value[2]
      self.__set__(self, instance, value[:2])
      return
    # checking for a 2-tuple to set envvar.
    if value != None and len(value) == 2:
      if hasattr(value[0], "__len__") and hasattr(value[1], "__len__"):
            is_2_tuple = len(value[0]) != 1 and len(value[1]) != 1
      else: is_2_tuple = value[0] == None or value[1] == None
      if is_2_tuple:
        self.envvar = value[1] 
        self.__set__(instance, value[0])
        return

    if value == None: self.directory = ""
    elif len(value.rstrip().lstrip()) == 0: self.directory = ""
    else: self.directory = relpath(value, self.fixed_point) 


