from abc import ABCMeta, abstractmethod
from collections import namedtuple

Program = namedtuple('Program', ['program', 'cmdline', 'directory', 'stdout', 'stderr'])
""" Holds data about program to launch. """

class Stop(StopIteration):
  """ No further calculation. """
  pass
class Fail(Exception):
  """ Process failed to run. """
  pass

class Process(object):
  """ Some methods and attributes process classes have in common. """
  __metaclass__ = ABCMeta
  def __init__(self, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    from copy import deepcopy
    from ..misc import RelativePath
    super(Process, self).__init__()
    self.jobfolder = jobfolder
    """ Folder to execute. """
    self.params = deepcopy(kwargs)
    """ Extra parameters passed on to functional's iterator. """
    self.nberrors = 0
    """ Number of restart on errors. """
    self.maxtrials = maxtrials
    """ Maximum number of restarts. """
    self.comm = comm
    """ MPI communicator. """

  @abstractmethod
  def poll(): 
    """ Polls current job. """
    pass

  def _cleanup(self):
    """ Cleanup files and crap. """
    if self.current_program[0] is None: return
    if not getattr(self.current_program[1].stdout, 'closed', True):
      self.current_program[1].stdout.close()
    if not getattr(self.current_program[1].stderr, 'closed', True): 
      self.current_program[1].stderr.close()
    if len(self.current_program) == 3 and self.current_program[2] is not None:
      from ..misc import LockFile
      LockFile(self.current_program[2].name).remove_stale()
      remove(self.current_program[2].name)
    # now remove reference to current program.
    self.current_program = tuple([None for i in self.current_program])
  def terminate(self):
    """ Terminates current process. """
    if self.current_program[0] is None: return
    try: self.current_program.terminate()
    except: pass
    self._cleanup()

  def kill(self):
    """ Kills current process. """
    if self.current_program[0] is None: return
    try: self.current_program.kill()
    except: pass
    self._cleanup()

  def _poll_current_process(self):
    """ Returns True if no error found.

        If no process is running, there is no error.
    """ 
    if self.current_program[0] is None: return True
    poll = self.current_program[0].poll()
    if poll is None: return
    self._cleanup()
    return poll < 0

  def wait(self):
    """ Waits for process to end, then cleanup. """
    if self.current_program[0] is None: return True
    return self.current_program[0].wait() >= 0



