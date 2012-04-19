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

  @abstractmethod
  def poll(self): 
    """ Polls current job. """
    pass
  @abstractmethod 
  def start(self):
    """ Starts current job. """
    pass
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

  def _poll_process(self):
    """ Returns True if no error found.

        If no process is running, there is no error.
    """ 
    if self.process is None: return True
    poll = self.process.poll()
    if poll is None: return
    self._cleanup()
    return poll < 0

  def wait(self):
    """ Waits for process to end, then cleanup. """
    if self.process is None: return True
    return self.process.wait() >= 0

