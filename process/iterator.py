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
    self._iterator = None
    """ Current iterator. """

  def poll(self):
    """ Polls current job. """
    from . import Fail

    # checks whether program was already started or not.
    if super(IteratorProcess, self).poll(): return True

    # check if we have currently running process.
    # catch StopIteration exception signaling that process finished.
    found_error = None
    try:
      if self.process.poll() == False: return False
    except Fail as failure: found_error = failure
    try: self.process._cleanup()
    finally: self.process = None

    # if stopped on or before previous job, then we have a retrial.
    if found_error is not None:
      self.nberrors += 1
      if self.nberrors >= self.maxtrials:
        self._cleanup()
        raise found_error if isinstance(found_error, Fail) \
              else Fail(str(found_error))

    # At this point, go to next iteration.
    process = self._get_process()
    if process is not None:
      self._next(process)
      return False

    self._cleanup()
    return True


  def _get_process(self):
    """ Finds next available iteration. """
    # first creates iterator, depending on input type.
    if self._iterator is None:
      iterator = self.functional.iter if hasattr(self.functional, 'iter')\
                 else self.functional
      self._iterator = iterator( comm=self._comm,
                                 outdir=self.outdir, 
                                 **self.params)
    try:
      result = self._iterator.next()
      while hasattr(result, 'success'): 
        result = self._iterator.next()
      return result
    except StopIteration: return None
    except Exception as e: raise Fail(e)

  def _next(self, process=None):
    """ Launches next process. """
    # start next process.
    self.process = process if process is not None else self._get_process()
    if self.process is None: return True
    self.process.start(self._comm.lend('all')) 
    return False

  def start(self, comm):
    """ Starts current job. """
    if super(IteratorProcess, self).start(comm): return True
    self._next()
    return False

  def wait(self):
    """ Waits for process to end, then cleanup. """
    if super(IteratorProcess, self).wait(): return True
    while not self.poll(): self.process.wait()
    self._cleanup()
    return False
