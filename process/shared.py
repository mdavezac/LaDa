from .process import Process
class SharedJobFolderProcess(IteratorProcess):
  """ Executes jobfolder with more than one job in child process. """
  def __init__(self, jobfolderpath, maxtrials=1, comm=None, nbpools=1, **kwargs):
    """ Initializes a process. """
    from copy import deepcopy
    super(JobFolderProcess, self).__init__(maxtrials, comm, **kwargs)
    self.nbpools = nbpools
    """ Number of pools over which to separate communicator. """
    self.jobfolderpath = jobfolderpath
    # start first process.
    self.poll()

  def poll(): 
    """ Polls current job. """
    from subprocess import Popen
    from shlex import split as split_cmd
    from misc import Changedir
    from . import Program
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    if self.current_program[0] is not None: 
      if self.current_program[0].poll() is None: return
      # kills stdout and stderr.
      if hasattr(self.current_program[1].stdout, 'close'): 
        self.current_program[1].stdout.close()
      if hasattr(self.current_program[1].stderr, 'close'): 
        self.current_program[1].stderr.close()
      # now remove reference to current program.
      self.current_program = None, None

    # At this point, loop until find something to do.
    found = False
    params = self.jobfolder.params.copy()
    params.update(self.params)
    for i, program in self.jobfolder.iterator(**params):
      if not getattr(program, 'success', False): 
        found = True
        break;
    # stop if no more jobs.
    if found == False: raise IteratorProcess.StopIteration()
    # if stopped on or before previous job, then we have a retrial.
    if i <= self.iter_index:
      if self.nberrors >= self.maxtrials: raise IteratorProcess.StopIteration()
      self.nberrors += 1
    # Open stdout and stderr if necessary.
    with Changedir(program.directory) as cwd:
     file_out = open(program.stdout, "a" if append else "w") \
                if program.stdout is not None else None 
     file_err = open(program.stderr, "a" if append else "w") \
                if program.stdout is not None else None 

    # now create process.
    program = Program(program.program, program.cmdline, program.directory, fileout, filerr)
    process = Popen(split_cmd(cmd), stdout=file_out, stderr=file_err, cwd=program.directory)
    self.current_program = process, program
