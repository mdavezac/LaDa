from .process import Process
class JobFolderProcess(Process):
  """ Executes folder in child process.
  
      Expects folder with a functional which does not have an iter method.
  """
  def __init__(self, jobfolder, outdir, maxtrials=1, comm=None, nbpools=1, **kwargs):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(JobFolderProcess, self).__init__(maxtrials, comm, **kwargs)

    self.jobfolder = jobfolder
    """ Job-folder to execute. """
    self.outdir = RelativePath(outdir).path
    """ Execution directory of the folder. """
    self.process = []
    """ List of currently running processes. """
    self.nbpools
    """ How many jobs to launch simultaneously. """
    self._finished = set()
    """ List of jobs to run. """
    self._torun = set(self.jobfolder.keys())
    """ List of errors. """
    self._errors = set()
    """ List of failed jobs. """

  @property
  def nbjobsleft(): 
    """ Number of jobs left. """
    return len(self._torun)

  def poll(): 
    """ Polls current job. """
    from tempfile import NamedTemporaryFile
    from os.path import join, exists, abspath
    from pickle import dump
    from sys import executable, pypath
    from ..misc import Program
    from .call import CallProcess
    from .iterator import IteratorProcess

    # check if we have currently running process.
    # catch StopIteration exception signaling that process finished.
    # catches Fail, but does nothing: each job is allowed to fail maxtrials
    # time.
    for i, (name, process) in enumerate(list(self.process)):
      try: delf.process.poll()
      except Stop: self.process.pop(i)
      except Fail as fail:
        self._errors.add(name)
        self.process.pop(i)

    # raises Stop if nothing left to do.
    if len(self.process) == 0 and len(self._torun) == 0:
      if len(self._errors) == 0: raise Stop()
      else: raise Fail()
    # Loop until all requisite number of processes is created, 
    # or until run out of jobs.
    while len(self.process) < self.nbpools and len(self._torun) > 0:
      name = self._torun.pop()
      # checks folder is still valid.
      if name not in self.jobfolder:
        from ..error import IndexError
        raise IndexError("Job-folder {0} no longuer exists.".format(name))
      jobfolder = self.jobfolder[name]
      if not jobfolder.is_executable:
        from ..error import IndexError
        raise IndexError("Job-folder {0} is no longuer executable.".format(name))
      # creates parameter dictionary. 
      params = self.jobfolder.params.copy()
      params.update(self.params)
      params['maxtrials'] = self.maxtrials
      params['comm'] = self.comm
      # chooses between an iterator process and a call process.
      if hasattr(jobfolder.functional, 'iter'):
        process = IteratorProcess(job.functional, join(self.outdir, name[1:]), **params)
      else:
        process = CallProcess(self.functional, join(self.outdir, name[1:]), **params)
      # appends process and starts it.
      self.process.append((name, process))
      process.start()

   def kill(self):
     """ Kills all currently running processes. """
     for name, process in self.process: process.kill()
   def terminate(self):
     """ Kills all currently running processes. """
     for name, process in self.process: process.terminate()
   def wait(self, sleep=0.5):
     """ Waits for all job-folders to execute and finish. """
     from time import sleep as time_sleep
     while True: 
       self.poll()
       time_sleep(sleep)
   def _cleanup(self):
     """ Cleans up after currently running processes. """
     for name, process in self.process: self.process._cleanup()
    
   def start(self):
     """ Start executing job-folders. """
     self.poll()
