""" Miscellaneous ressources. """
__all__ = ['Program', 'execute_program', 'copyfile', 'changedir', \
           'read_input', 'exec_input', 'load', 'RelativePath']
from collections import namedtuple
from types import ModuleType

from changedir import Changedir
from relativepath import RelativePath

Program = namedtuple('Program', ['program', 'cmdline', 'directory', 'stdout', 'stderr'])
""" Holds data about program to launch. """

def execute_program(program, append=False, mpirun_exe=None, comm=None, **kwargs):
   """ Executes a program.
   
       Launches a program, generally using the :py:data:`mpirun_exe
       <lada.mpirun_exe>` format string.

       :param program:
          Contains the data relevant to launching an external program.
       :type program: :py:Class:`Program`
       :param Boolean append:
          Whether to create new standard output and error, or whether to
          append.
       :param mpirun_exe:
          String to use as format for the full commandline. If None, uses
          :py:data:`mpirun_exe <lada.mpirun_exe>` by default.
       :param comm:
          Dictionary containing information about launching an mpi program.
          If None, defaults to :py:data:`default_comm <lada.default_comm>`.
   """
   from subprocess import Popen
   from shlex import split as split_cmd
   from lada import mpirun_exe as mpirun_exe_global, default_comm
   
   d = {}
   d['program'] = program.program
   d['cmdline'] = program.cmdline 
   d.update(comm if comm is not None else default_comm)
   d.update(kwargs)
   if mpirun_exe == None: mpirun_exe = mpirun_exe_global
   cmd = mpirun_exe.format(**d)
   
   with Changedir(program.directory) as cwd:
     file_out = open(program.stdout, "a" if append else "w") \
                if program.stdout is not None else None 
     file_err = open(program.stderr, "a" if append else "w") \
                if program.stdout is not None else None 
     try:
       vasp_proc = Popen(split_cmd(cmd), stdout=file_out, stderr=file_err)
       vasp_proc.wait()
     finally:
       file_out.close()
       file_err.close()


def copyfile(src, dest=None, nothrow=None, symlink=False, aslink=False, nocopyempty=False):
  """ Copy ``src`` file onto ``dest`` directory or file.

      :param src:
          Source file.
      :param dest: 
          Destination file or directory.
      :param nothrow:
          Throwing is disable selectively depending on the content of nothrow:
          - *exists*: will not throw is src does not exist.
          - *isfile*: will not throw is src is not a file.
          - *same*: will not throw if src and dest are the same.
          - *none*: ``src`` can be None.
          - *null*: ``src`` can be '/dev/null'.
          - *never*: will never throw.

      :param symlink:
          Creates link rather than actual hard-copy. Symlink are
          created with relative paths given starting from the directory of
          ``dest``.  Defaults to False.
      :param aslink: 
          Creates link rather than actual hard-copy *if* ``src`` is
          itself a link. Links to the file which ``src`` points to, not to
          ``src`` itself. Defaults to False.
      :parma nocopyempty:
          Does not perform copy if file is empty. Defaults to False.

      This function fails selectively, depending on what is in ``nothrow`` list.
  """
  try:
    from os import getcwd, symlink as ln, remove
    from os.path import isdir, isfile, samefile, exists, basename, dirname,\
                        join, islink, realpath, relpath, getsize
    from shutil import copyfile as cpf
    # sets up nothrow options.
    if nothrow is None: nothrow = []
    if isinstance(nothrow, str): nothrow = nothrow.split()
    if nothrow == 'all': nothrow = 'exists', 'same', 'isfile', 'none', 'null'
    nothrow = [u.lower() for u in nothrow]
    # checks and normalizes input.
    if src is None: 
      if 'none' in nothrow: return False
      raise IOError("Source is None.")
    if dest is None: dest = getcwd()
    if dest == '/dev/null': return True
    if src  == '/dev/null':
      if 'null' in nothrow: return False
      raise IOError("Source is '/dev/null' but Destination is {0}.".format(destination))

    # checks that input source file exists.
    if not exists(src): 
      if 'exists' in nothrow: return False
      raise IOError("{0} does not exist.".format(src))
    if not isfile(src):
      if 'isfile' in nothrow: return False
      raise IOError("{0} is not a file.".format(src))
    # makes destination a file.
    if exists(dest) and isdir(dest): dest = join(dest, basename(src))
    # checks if destination file and source file are the same.
    if exists(dest) and samefile(src, dest): 
      if 'same' in nothrow: return False
      raise IOError("{0} and {1} are the same file.".format(src, dest))
    if nocopyempty and isfile(src):
      if getsize(src) == 0: return
    if aslink and islink(src): symlink, src = True, realpath(src)
    if symlink:
      if exists(dest): remove(dest)
      if relpath(src, dirname(dest)).count("../") == relpath(src, '/').count("../"):
        ln(src, realpath(dest))
      else:
        with Changedir(dirname(dest)) as cwd:
           ln(relpath(src, dirname(dest)), basename(dest))
    else: cpf(src, dest)
  except:
    if 'never' in nothrow: return False
    raise
  else: return True

class Input(ModuleType):
  """ Fake class which will be updated with the local dictionary. """
  def __init__(self, name = "lada_input"): 
    """ Initializes input module. """
    super(Input, self).__init__(name, "Input module for lada scripts.")
  def __getattr__(self, name):
    raise AttributeError( "All out of cheese!\n"
                          "Required input parameter '{0}' not found in {1}." \
                          .format(name, self.__name__) )
  def __delattr__(self, name):
    raise RuntimeError("Cannot delete object from input namespace.")
  def __setattr__(self, name, value):
    raise RuntimeError("Cannot set/change object in input namespace.")
  def update(self, other):
    if hasattr(other, '__dict__'): other = other.__dict__
    for key, value in other.items():
      if key[0] == '_': continue
      super(Input, self).__setattr__(key, value)
  @property
  def __all__(self):
    return list([u for u in self.__dict__.iterkeys() if u[0] != '_'])
  def __contains__(self, name):
    return name in self.__dict__

def read_input(filename='input.py', global_dict=None, local_dict = None, paths=None, comm = None):
  """ Reads and executes input script and returns local dictionary (as namespace instance). """
  from os.path import exists, basename
  assert exists(filename), IOError('File {0} does not exist.'.format(filename))
  with open(filename, 'r') as file: string = file.read()
  return exec_input(string, global_dict, local_dict, paths, comm, basename(filename))

def exec_input( script, global_dict=None, local_dict = None,\
                paths=None, comm = None, name=None):
  """ Executes input script and returns local dictionary (as namespace instance). """
  # stuff to import into script.
  from os import environ
  from os.path import abspath, expanduser
  from math import pi 
  from numpy import array, matrix, dot, sqrt, abs, ceil
  from numpy.linalg import norm, det
  from .. import physics, crystal
  from . import Input
  
  # Add some names to execution environment.
  if global_dict is None: global_dict = {}
  global_dict.update( { "environ": environ, "pi": pi, "array": array, "matrix": matrix, "dot": dot,
                        "norm": norm, "sqrt": sqrt, "ceil": ceil, "abs": abs,  "det": det,
                        "physics": physics, "expanduser": expanduser, "load": load })
  for key in crystal.__all__: global_dict[key] = getattr(crystal, key)
  if local_dict is None: local_dict = {}
  # Executes input script.
  exec(script, global_dict, local_dict)

  # Makes sure expected paths are absolute.
  if paths is not None:
    for path in paths:
      if path not in local_dict: continue
      local_dict[path] = abspath(expanduser(local_dict[path]))
    
  if name is None: name = 'None'
  result = Input(name)
  result.update(local_dict)
  return result

def load(data, *args, **kwargs):
  """ Loads data from the data files. """
  from os import environ
  from os.path import dirname, exists, join
  if "directory" in kwargs: 
    raise KeyError("directory is a reserved keyword of load")

  # find all possible data directories
  directories = []
  if "data_directory" in globals():
    directory = globals()["data_directory"]
    if hasattr(directory, "__iter__"): directories.extend(directory)
    else: directories.append(directory)
  if "LADA_DATA_DIRECTORY" in environ:
    directories.extend(environ["LADA_DATA_DIRECTORY"].split(":"))

  # then looks for data file.
  if data.rfind(".py") == -1: data += ".py"
  for directory in directories:
    if exists(join(directory, data)):
      kwargs["directory"] = dirname(join(directory, data))
      result = {}
      execfile(join(directory, data), {}, result)
      return result["init"](*args, **kwargs)
  raise IOError("Could not find data ({0}).".format(data))

def add_setter(method, docstring = None): 
  """ Adds an input-like setter property. """
  def _not_available(self): raise RuntimeError("Error: No cheese available.")
  if docstring is None and hasattr(method, "__doc__"): docstring = method.__doc__
  return property(fget=_not_available, fset=method,  doc=docstring)
