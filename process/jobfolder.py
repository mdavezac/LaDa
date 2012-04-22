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
    self.nbpools = nbpools
    """ How many jobs to launch simultaneously. """
    self._finished = set()
    """ List of jobs to run. """
    self._torun = set(self.jobfolder.keys())
    """ List of errors. """
    self.errors = {}
    """ Map between name of failed jobs and exception. """

  @property
  def nbjobsleft(self): 
    """ Number of jobs left. """
    return len(self._torun)

  def poll(self): 
    """ Polls current job. """
    from os.path import join
    from ..error import IndexError
    from . import Fail
    from .call import CallProcess
    from .iterator import IteratorProcess

    if len(self.process) == 0  and self.started: return True
    self.started = True
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
      self.process.pop(i)[1]._cleanup()

    # returns True if nothing left to do.
    if len(self.process) == 0 and len(self._torun) == 0:
      if len(self.errors) == 0: return True
      else: raise Fail()
    # Loop until all requisite number of processes is created, 
    # or until run out of jobs.
    while len(self.process) < self.nbpools and len(self._torun) > 0:
      name = self._torun.pop()
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
      params['comm'] = self.comm
      # chooses between an iterator process and a call process.
      if hasattr(jobfolder.functional, 'iter'):
        process = IteratorProcess(jobfolder.functional, join(self.outdir, name), **params)
      else:
        process = CallProcess(self.functional, join(self.outdir, name), **params)
      # appends process and starts it.
      self.process.append((name, process))
      process.start()
    return False

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
  def wait(self, sleep=0.5):
    """ Waits for all job-folders to execute and finish. """
    from time import sleep as time_sleep
    while not self.poll(): 
      if sleep > 0: time_sleep(sleep)
  def _cleanup(self):
    """ Cleans up after currently running processes. """
    for name, process in self.process: self.process._cleanup()
   
  def start(self):
    """ Start executing job-folders. """
    self.poll()
