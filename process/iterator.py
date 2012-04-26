from .process import Process
class IteratorProcess(Process):
  """ Executes an iteration function in child process. 
  
      Expects folder with a functional which has an iter method.
      This method should yield either an extractor object for previously
      successful runs, or a derived :py:class:`~lada.process.process.Process`
      instance to execute.
  """
  def __init__(self, functional, outdir, maxtrials=1, **kwargs):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(IteratorProcess, self).__init__(maxtrials, **kwargs)
    self.functional = functional
    """ Iterable to execute. """
    self.outdir = RelativePath(outdir).path
    """ Execution directory of the folder. """
    self._iterindex = 0
    """ Current iteration. """

  def poll(self, wait=False): 
    """ Polls current job. """
    from . import Fail

    # checks whether program was already started or not.
    if super(IteratorProcess, self).poll(): return True

    # check if we have currently running process.
    # catch StopIteration exception signaling that process finished.
    found_error = None
    if self.process is not None:
      try:
        if wait == True: self.process.wait()
        elif self.process.poll() == False: return False
      except Fail as failure: found_error = failure
      self._cleanup()
    # At this point, loop until find something to do.
    found = False
    if hasattr(self.functional, 'iter'):
      def iterator(**params):
        for dummy in self.functional.iter(**params): yield dummy
    else: iterator = self.functional
    for i, process in enumerate(iterator(comm=self._comm, outdir=self.outdir, **self.params)):
      if not getattr(process, 'success', False): 
        found = True
        break;
    # stop if no more jobs.
    if found == False: return True
    # if stopped on or before previous job, then we have a retrial.

    if i <= self._iterindex and found_error is not None:
      self.nberrors += 1
      if self.nberrors >= self.maxtrials:
        raise found_error if isinstance(found_error, Fail) \
              else Fail(str(found_error))
    self._iterindex = i

    # start next process.
    self.process = process
    self.process.start() 
    return False

  def start(self, comm):
    """ Starts current job. """
    super(IteratorProcess, self).start(comm)
    self.poll()

  def wait(self):
    """ Waits for process to end, then cleanup. """
    if self.process is None:
      if self.started: return True
      self.start()
    return self.poll(wait=True)
