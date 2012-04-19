class IteratorProcess(Process):
  """ Executes an iteration function in child process. 
  
      Expects folder with a functional which has an iter method.
      This method should yield either an extractor object for previously
      successful runs, or a :py:class:`misc.Program` object to execute.
  """
  def __init__(self, functional, outdir, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(IteratorProcess, self).__init__(maxtrials, comm, **kwargs)
    self.functional = functional
    """ Folder to execute. """
    self.outdir = RelativePath(outdir).path
    """ Execution directory of the folder. """
    self._iterindex = 0
    """ Current iteration. """
    self.poll()

  def poll(): 
    """ Polls current job. """
    from . import Program
    from .process import Stop, Fail
    # check if we have currently running process.
    # catch StopIteration exception signaling that process finished.
    found_error = None
    if self.process is not None:
      try: self.process.poll()
      except Stop: pass
      except Fail as fail: found_error = fail
    # At this point, loop until find something to do.
    found = False
    for i, process in self.functional(self.comm=comm, **self.params):
      if not getattr(process, 'success', False): 
        found = True
        break;
    # stop if no more jobs.
    if found == False: raise Stop()
    # if stopped on or before previous job, then we have a retrial.

    if i <= self._iterindex and found_error is not None:
      self.nberrors += 1
      if self.nberrors >= self.maxtrials: raise Fail()
    self._iterindex = i

    # start next process.
    self.process = process
    self.process.start() 

  def start(self):
    """ Starts current job. """
    self.poll()
