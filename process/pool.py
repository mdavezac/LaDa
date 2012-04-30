from .jobfolder import JobFolderProcess
class PoolProcess(JobFolderProcess):
  """ Executes folder in child processes.
  
      Expects a jobfolder on input. Executable job-folders are launched in
      parallel, with up to :py:attr:`~PoolProcess.nbpools` running
      instances. Each instance is allocated an equal number of machines.
  """
  def __init__(self, jobfolder, outdir, maxtrials=1, processalloc=None, **kwargs):
    """ Initializes a process. """
    super(PoolProcess, self).__init__(jobfolder, outdir, maxtrials, **kwargs)
    del self.nbpools # not needed here.

    self.processalloc = processalloc
    """ How many jobs to launch simultaneously. """

  @property
  def nbjobsleft(self): 
    """ Number of jobs left. """
    return len(self._torun)

  def poll(self): 
    """ Polls current job. """
    from . import Fail

    if super(PoolProcess, self).poll(): return True

    # weed out successfull and failed jobs.
    finished = []
    for i, (name, process) in enumerate(list(self.process)):
      try:
        if process.poll() == True: 
          self._finished.add(name)
          finished.append(i)
      except Fail as fail:
        self.errors[name] = fail
        finished.append(i)
    for i in sorted(finished)[::-1]:
      name, process = self.process.pop(i)
      process._cleanup()

    # returns True if nothing left to do.
    if len(self.process) == 0 and len(self._torun) == 0:
      if len(self.errors) == 0: return True
      else: raise Fail()

    # adds new jobs. 
    self._next()
    return False
 
  def _next(self):
    """ Adds more processes.
    
        This is the subroutine to overload in a derived class which would
        implement some sort of scheduling.
    """
    from os.path import join
    from ..error import IndexError
    from .call import CallProcess
    from .iterator import IteratorProcess

    # nothing else to do.
    if len(self._torun) == 0: return
    # cannot add more processes.
    if len(self.process) >= self.nbpools: return
    # no more machines to allocate...
    if self._comm['n'] == 0: return

    # split processes into local comms. Make sure we don't oversuscribe.
    jobs = self._getjobs()
    try: 
      # Loop until all requisite number of processes is created, 
      # or until run out of jobs, or until run out of comms. 
      for name in jobs:
        self._torun.pop(self._torun.index(name))
        nprocs = self._alloc.pop(name)
        # checks folder is still valid.
        if name not in self.jobfolder:
          raise IndexError("Job-folder {0} no longuer exists.".format(name))
        jobfolder = self.jobfolder[name]
        if not jobfolder.is_executable:
          raise IndexError("Job-folder {0} is no longuer executable.".format(name))
        # creates parameter dictionary. 
        params = jobfolder.params.copy()
        params.update(self.params)
        params['maxtrials'] = self.maxtrials
        # chooses between an iterator process and a call process.
        if hasattr(jobfolder.functional, 'iter'):
          process = IteratorProcess(jobfolder.functional, join(self.outdir, name), **params)
        else:
          process = CallProcess(self.functional, join(self.outdir, name), **params)
        # appends process and starts it.
        self.process.append((name, process))
        process.start(self._comm.acquire(nprocs))
    except:
      self.terminate()
      raise

  def kill(self):
    """ Kills all currently running processes. """
    for name, process in self.process: process.kill()
    self._cleanup()
    self.process = []
    
  def terminate(self):
    """ Kills all currently running processes. """
    for name, process in self.process: process.terminate()
    self._cleanup()
    self.process = []

  @property 
  def done(self):
    """ True if job already finished. """
    return self.started and len(self.process) == 0
  @property
  def nbrunning_processes(self):
    """ Number of running processes. 

        For simple processes, this will be one or zero.
        For multitasking processes this may be something more.
    """
    return 0 if (not self.started) or len(self.process) == 0 else 1

  def wait(self, sleep=0.5):
    """ Waits for all job-folders to execute and finish. """
    from time import sleep as time_sleep
    if super(PoolProcess, self).wait(): return True
    while not self.poll(): 
      if sleep > 0: time_sleep(sleep)
    return False

  def _cleanup(self):
    """ Cleans up after currently running processes. """
    try: 
      for name, process in self.process: process._cleanup()
    finally:
      self.process = []
      if hasattr(self, '_comm'): 
        try: self._comm.cleanup()
        finally: del self._comm

  def _getjobs(self):
    """ List of jobs to run. """
    from operator import itemgetter
    N = self._comm['N']

    # creates list of possible jobs.
    availables = sorted( [(key, u) for key, u in self._alloc.iteritems() if u <= N],
                         key=itemgetter(1) )

    # hecks first if any jobs fits exactly the available number of nodes.
    if availables[-1][1] == N: return [availables[-1][1]]

    # creates a map of bins: 
    bins = {}
    for key, u in availables: 
      if u not in bins: bins[u] = 1
      else: bins[u] += 1

    def get(n, bins, xvec):
      """ Loops over possible combinations. """
      from random import choice
      key = choice(list(bins.keys()))
      for u in xrange(min(bins[key], n // key), -1, -1):
        newbins = bins.copy()
        del newbins[key]
        newn = n - u * key
        if newn == 0:
          yield xvec + [(key, u)], True
          break
        for key in list(newbins.keys()):
          if newbins[key] > newn: del newbins[key]
        if len(newbins) == 0:
          yield xvec + [(key, u)], False
          continue
        for othervec, perfect in get(newn, newbins, xvec + [(key, u)]):
          yield othervec, perfect
          if perfect: return

    xvec = []
    nprocs, njobs = 0, 0
    for u, perfect in get(N, bins, []):
      if perfect: xvec = u; break
      p, j = sum(a*b for a, b in u), sum(a for a, b in u)
      if p > nprocs or (p == nprocs and j < njobs):
        xvec, nprocs, njobs = list(u), p, j

    # now we have a vector with the right number of jobs, but not what those
    # jobs are.
    results = []
    for key, value in xvec:
      withkeyprocs = [name for name, n in availables if n == key]
      results.extend(withkeyprocs[:value])
    return results
    
   
  def start(self, comm):
    """ Start executing job-folders. """
    from .process import Process
    from .mpi import MPISizeError
    if isinstance(self.processalloc, int): 
      self.nbpools = comm['n'] // self.processalloc
      return super(PoolProcess, self).start(comm)

    if Process.start(self, comm): return True

    # maps jobs and requested allocation
    self._alloc = {}
    for key, job in self.jobfolder.iteritems(): 
      self._alloc[key] = self.mppalloc(job)

    # check max job size.
    toolarge = [key for key, u in self._alloc.iteritems() if u > comm['u']]
    if len(toolarge):
      raise MPISizeError( "The following jobs require too many processors:\n"\
                                "{0}\n".format(toolarge) )

    self._next()
    return False

