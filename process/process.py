from abc import ABCMeta, abstractmethod

class Process(object):
  """ Some methods and attributes process classes have in common. """
  __metaclass__ = ABCMeta
  def __init__(self, maxtrials=1, **kwargs):
    """ Initializes a process. """
    super(Process, self).__init__()

    self.params = kwargs
    """ Extra parameters passed on to functional's iterator. """
    self.nberrors = 0
    """ Number of restart on errors. """
    self.maxtrials = maxtrials
    """ Maximum number of restarts. """
    self.process = None
    """ Currently running process. """
    self.started = False
    """ Whether the program was ever started. """

  @abstractmethod
  def poll(self): 
    """ Polls current job. """
    from ..error import internal
    if not self.started: raise internal("Process was never started.")
    return self.nbrunning_processes == 0

  @property 
  def done(self):
    """ True if job already finished. """
    return self.started and self.process is None

  @property
  def nbrunning_processes(self):
    """ Number of running processes. 

        For simple processes, this will be one or zero.
        For multitasking processes this may be something more.
    """
    return 0 if (not self.started) or self.process is None else 1

  @abstractmethod
  def start(self, comm):
    """ Starts current job. """
    from . import AlreadyStarted
    from .mpi import MPISizeError, Communicator
    if self.done: return True
    if self.started: raise AlreadyStarted('start cannot be called twice.')
    self.started = True
    if comm is not None:
      if comm['n'] == 0: raise MPISizeError('Empty communicator passed to process.')
      self._comm = comm if hasattr(comm, 'machines') else Communicator(**comm) 
    return False

  def _cleanup(self):
    """ Cleans up behind process instance.
    
        This may mean closing standard input/output file, removing temporary
        files.. By default, calls cleanup of process, and sets process to None.
    """
    try:
      if hasattr(self.process, '_cleanup'): self.process._cleanup()
    finally:
      self.process = None
      if hasattr(self, '_comm'): 
        try: self._comm.cleanup()
        finally: del self._comm

  def terminate(self):
    """ Terminates current process. """
    if self.process is None: return
    try: self.process.terminate()
    except: pass
    self._cleanup()

  def kill(self):
    """ Kills current process. """
    if self.process is None: return
    try: self.process.kill()
    except: pass
    self._cleanup()

  @abstractmethod
  def wait(self):
    """ Waits for process to end, then cleanup. """
    from ..error import internal
    if not self.started: raise internal("Process was never started.")
    if self.nbrunning_processes == 0: return True
