""" Miscellaneous """
import _opt 
from contextlib import contextmanager
__load_vasp_in_global_namespace__ = _opt.__load_vasp_in_global_namespace__
__load_pescan_in_global_namespace__ = _opt.__load_pescan_in_global_namespace__
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
    except: print "Am here"
    if self.unit == streams.input:    os.dup2(self.old, sys.__stdin__.fileno()) 
    elif self.unit == streams.output: os.dup2(self.old, sys.__stdout__.fileno())
    elif self.unit == streams.error:  os.dup2(self.old, sys.__stderr__.fileno())
    else: raise RuntimeError("Unknown redirection unit.")

def redirect_all(output=None, error=None, input=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param output: Filename to which to redirect output. 
      @param error: Filename to which to redirect err. 
      @param input: Filename to which to redirect input. 
      @param append: If true, will append to files. All or nothing.
  """
  from contextlib import nested
  result = []
  for value, unit in [ (output, streams.output), (error, streams.error), (input, streams.input) ]:
    if value == None: continue
    result.append( _RedirectAll(unit=unit, filename=value, append=append) )
  return nested(*result)


def redirect(fout=None, ferr=None, fin=None, pyout=None, pyerr=None, pyin=None, append = False):
  """ A context manager to redirect inputs, outputs, and errors. 
  
      @param fout: Filename to which to redirect fortran output. 
      @param ferr: Filename to which to redirect fortran err. 
      @param fin: Filename to which to redirect fortran input. 
      @param pyout: Filename to which to redirect C/C++ output. 
      @param pyerr: Filename to which to redirect C/C++ err. 
      @param pyin: Filename to which to redirect C/C++ input. 
      @param append: If true, will append to files. All or nothing.
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
  class Input:
    def __init__(self): self._inputfilepath = filename
    def __getattr__(self, name):
      raise AttributeError( "All out of cheese!\nRequired input parameter \"%s\" not found in %s." \
                            % (name, self._inputfilepath) )
  result = Input()
  result.__dict__.update(local_dict)
  return result


  
class LockFile(object):
  """ Gets an advisory lock for file C{filename}.

      Creates a lock directory named after a file. Relies on the presumably
      atomic nature of creating a directory. 
  """
  def __init__(self, filename, timeout = None, sleep = 0.5):
    """ Creates a lock object. 

        Does not acquire lock at this stage. 
        @param timeout: will raise a RuntimeError when calling L{self.lock} if
          the lock could not be aquired within this time.
        @param sleep: Time to sleep between checks when waiting to acquire lock. 
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
    from os import makedirs, error
    import time

    start_time = time.time()
    # loops until acqires lock.
    while self._owns_lock == False: 
      # loops until directory ceases to exist.
      while self.is_locked: 
        time.sleep(self.sleep)
        if self.timeout != None:
          assert time.time() - start_time < self.timeout, \
                 RuntimeError("Could not acquire lock on %s." % (filename))
      # tries to create director.
      try: makedirs(self.lock_directory)
      # if fails, then another process already created it. Just keep looping.
      except error: continue
      # Done!
      else: self._owns_lock = True


  @property
  def lock_directory(self):
    """ Name of lock directory. """
    from os.path import abspath, dirname, join, basename
    return join(dirname(abspath(self.filename)), "." + basename(self.filename) + "-lada_lockdir")

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

@contextmanager
def acquire_lock(filename, timeout = None, sleep = 0.5):
  """ Context for a L{LockFile}. """
  # Enter context.
  lock = LockFile(filename, timeout, sleep)
  lock.lock()
  yield lock
  # exit context.
  lock.release()

@contextmanager
def open_exclusive(filename, mode="r", timeout = None, sleep = 0.5):
  """ Opens file while checking for advisory lock.

      This context uses fcntl to acquire a lock on file before reading or
      writing. This is an advisory lock only.
  """
  # Enter context.
  lock = LockFile(filename, timeout, sleep)
  lock.lock()
  file = open(filename, mode)
  yield file
  # Exit context.
  file.close()
  lock.release()

