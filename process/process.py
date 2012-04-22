from abc import ABCMeta, abstractmethod

class Process(object):
  """ Some methods and attributes process classes have in common. """
  __metaclass__ = ABCMeta
  def __init__(self, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    super(Process, self).__init__()

    self.params = kwargs
    """ Extra parameters passed on to functional's iterator. """
    self.nberrors = 0
    """ Number of restart on errors. """
    self.maxtrials = maxtrials
    """ Maximum number of restarts. """
    self.comm = comm
    """ MPI communicator. """
    self.process = None
    """ Currently running process. """
    self.started = False
    """ Whether the program was ever started. """

  @abstractmethod
  def poll(self): 
    """ Polls current job. """
    if self.started and self.process is None: return True
    self.started = True

  @abstractmethod
  def start(self):
    """ Starts current job. """
    if self.started and self.process is None: return True
    self.started = True
    return False

  def _cleanup(self):
    """ Cleans up behind process instance.
    
        This may mean closing standard input/output file, removing temporary
        files.. By default, calls cleanup of process, and sets process to None.
    """
    try:
      if hasattr(self.process, '_cleanup'): self.process._cleanup()
    finally: self.process = None

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
    pass
