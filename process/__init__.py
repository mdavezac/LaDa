""" Manages external programs. """
__docformat__ = "restructuredtext en"
__all__  = [ 'Process', 'ProgramProcess', 'CallProcess', 'IteratorProcess',
             'JobFolderProcess', 'PoolProcess', 'Fail', 'which', 'DummyProcess' ]

from ..error import root
from pool import PoolProcess
from process import Process
from call import CallProcess
from program import ProgramProcess
from iterator import IteratorProcess
from jobfolder import JobFolderProcess
from dummy import DummyProcess

class ProcessError(root):
  """ Root of special exceptions issued by process module. """
class Fail(ProcessError):
  """ Process failed to run successfully. """
  pass
class AlreadyStarted(ProcessError):
  """ Process already started. """

def which(program):
  """ Gets location of program by mimicking bash which command. """
  from os import environ, getcwd
  from os.path import split, expanduser, expandvars, join
  from itertools import chain
  from ..misc import RelativePath
  from ..error import IOError

  def is_exe(path):
    from os import access, X_OK
    from os.path import isfile
    return isfile(path) and access(path, X_OK)

  exprog = expanduser(expandvars(program))
  fpath, fname = split(exprog)
  if fpath:
    if is_exe(exprog): return RelativePath(exprog).path
  else:
    for dir in chain([getcwd()], environ["PATH"].split(':')):
      if is_exe(join(dir, exprog)): return RelativePath(join(dir, exprog)).path

  raise IOError('Could not find executable {0}.'.format(program))

