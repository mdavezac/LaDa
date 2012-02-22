""" Miscellaneous ressources. """
__all__ = ['Program', 'execute_program', 'add_setter', 'copyfile']
from collections import namedtuple

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

