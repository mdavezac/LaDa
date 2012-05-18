from .jobfolder import JobFolderProcess
class PoolProcess(JobFolderProcess):
  """ Executes folder in child processes.
  
      Expects a jobfolder on input. Executable job-folders are launched in
      parallel, with up to :py:attr:`~PoolProcess.nbpools` running
      instances. Each instance is allocated an equal number of machines.
  """
  def __init__(self, jobfolder, outdir, processalloc, maxtrials=1, **kwargs):
    """ Initializes a process. """
    super(PoolProcess, self).__init__(jobfolder, outdir, maxtrials, **kwargs)
    del self.nbpools # not needed here.

    self.processalloc = processalloc
    """ How many jobs to launch simultaneously. """
    self._alloc = {}
    """ Maps job vs rquested process allocation. """
    for name in self._torun:
      self._alloc[name] = self.processalloc(self.jobfolder[name])
    assert len(set(self._alloc.keys())) == len(self._alloc)


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
    # no more machines to allocate...
    if self._comm['n'] == 0: return
    # cannot add more processes.
    if isinstance(self.processalloc, int):
      return super(PoolProcess, self)._next()

    # split processes into local comms. Make sure we don't oversuscribe.
    jobs = self._getjobs()
    assert sum(self._alloc[u] for u in jobs) <= self._comm['n']
    try: 
      # Loop until all requisite number of processes is created, 
      # or until run out of jobs, or until run out of comms. 
      for name in jobs:
        self._torun = self._torun - set([name])
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
        try: process.start(self._comm.lend(nprocs))
        except Exception as e:
          self.errors[name] = e
          name, process = self.process.pop(i)
          process._cleanup()
    except:
      self.terminate()
      raise

  def _getjobs(self):
    """ List of jobs to run. """
    from operator import itemgetter
    N = self._comm['n']

    # creates list of possible jobs.
    availables = sorted( [(key, u) for key, u in self._alloc.iteritems() if u <= N],
                         key=itemgetter(1) )
    if len(availables) == 0: return []

    # hecks first if any jobs fits exactly the available number of nodes.
    if availables[-1][1] == N: return [availables[-1][0]]

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
        for v in list(newbins.keys()):
          if v > newn: del newbins[v]
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

    # check max job size.
    toolarge = [key for key, u in self._alloc.iteritems() if u > comm['n']]
    if len(toolarge):
      raise MPISizeError( "The following jobs require too many processors:\n"\
                                "{0}\n".format(toolarge) )

    self._next()
    return False

  def update(self, jobfolder, deleteold=False):
    """ Updates list of jobs.
    
        Adds jobfolders which are not in ``self.jobfolder`` but in the input.
        Deletes those which in ``self.jobfolder`` but not in the input.
        Does nothing if job is currently running.
        Finished jobs are not updated.
        If ``deleteold`` is True, then removed finished jobs from job-folder.

        Processes jobfolder from root, even if passed a child folder.
    """
    running = set([n for n in self.process])
    for name, value in jobfolder.root.iteritems():
      if name in running: continue
      elif name not in self.jobfolder.root:
        newjob = self.jobfolder.root / name
        newjob.functional = value.functional
        newjob.params.update(value.params)
        for key, value in value.__dict__.iteritems():
          if key in ['children', 'params', '_functional', 'parent']: continue
          setattr(self, key, value)
        self._torun.add(name)
        self._alloc[name] = self.processalloc(self.jobfolder.root[name])
      elif name not in self._finished:
        self.jobfolder.root[name] = value
        self._alloc[name] = self.processalloc(self.jobfolder.root[name])
    for name in self.jobfolder.root.iterkeys():
      if name in self._finished and deleteold:
        del self.jobfolder.root[name]
        self._alloc.pop(name)
      elif name not in jobfolder.root:
        if name in running: continue
        del self.jobfolder.root[name]
        if name in self._torun: self._torun.remove(name)
        self._alloc.pop(name, None)



