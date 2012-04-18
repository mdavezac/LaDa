""" Functor to *bleed* a job folder amongst pools of processes.
    
    Bleeding a job folder means that pools of processes are actively
    picking jobs from the dictionary, marking them as executed, and running
    them. Any one job should only be visited once.
"""
__docformat__ = "restructuredtext en"
__all__ = ['findone', 'iter', 'execall']
from abc import ABCMeta, abstractmethod


@contextmanager
def find_one(path, jobname=None):
  """ Returns an executabe job-folder as a context. 

      The goal is to retrieve an executable job-folder from disk, while using a
      lock-file to prevent other processes from reading/writing the job-folder
      at an importune time. The lock-file is not held throughout the existence
      of the context. It is only held for specific operations when entering and
      leaving the context. The function is a context_. As such, it is possible
      to modify the job-folder. The modifications will be saved when leaving
      the context.

      .. _context: http://docs.python.org/reference/datamodel.html#context-managers
  """ 
  from ..misc import LockFile

  found = False
  # acquire lock. Yield context outside of lock!
  with LockFile(path) as lock:
    # Loads pickle.
    with open(self._filename, 'r') as file: jobfolder = pickle_load(file)
    # Finds first untagged job.
    if jobname is None:
      for job in jobfolder.itervalues():
        if not job.is_tagged:
          found = True
          break
    else: job = jobfolder[jobname]
    # tag the job before surrendering the lock.
    if found and jobname is None:
      job.tag()
      with open(self._filename, 'w') as file: dump(jobfolder, file)

  # Check we found an untagged job. Otherwise, we are done.
  if not found:
    yield None
    return
    
  # context returns jobs.
  yield job

  # saves job since it might have been modified
  # acquire a lock first.
  with LockFile(self._filename) as lock:
    # Loads pickle.
    with open(self._filename, 'r') as file: jobfolder = pickle_load(file)
    # modifies job.
    jobfolder[job.name] = job
    # save jobs.
    with open(self._filename, 'w') as file: dump(jobfolder, file)

def iterall(path):
  """ Iterates over executable job-folders on disk. 
  
      The job-folder is continuously read from disk (and locked when
      read/written). This way, many processes can share access to the same
      job-folder.
  """
  from tempfile import NamedTemporaryFile
  from pickle import dump
  from ..opt import RelativeDirectory

  from os.path import exists
  from pickle import load as pickle_load, dump
  from ..misc import LockFile
  from ..error import IOError

  # tries to read file.
  if not exists(path): raise IOError('Job-folder does not exist.')
  # infinite loop. breaks when no new jobs can be found.
  while True:
    # only local root reads stuff. 
    job = None
    with findone(path) as job:
      if job is None: break
      yield job

class Process(object):
  """ Some methods and attributes processes have in common. """
  __metaclass__ = ABCMeta
  class StopIteration(StopIteration):
    """ No further iteration. """
    pass
  class Failed(Exception):
    """ Process failed to run. """
    pass

  def __init__(self, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    from copy import deepcopy
    from ..misc import RelativePath
    super(Process, self).__init__()
    self.jobfolder = jobfolder
    """ Folder to execute. """
    self.params = deepcopy(kwargs)
    """ Extra parameters passed on to functional's iterator. """
    self.nberrors = 0
    """ Number of restart on errors. """
    self.maxtrials = maxtrials
    """ Maximum number of restarts. """
    self.comm = comm
    """ MPI communicator. """

  @abstractmethod
  def poll(): 
    """ Polls current job. """
    pass

  def _cleanup(self):
    """ Cleanup files and crap. """
    if self.current_program[0] is None: return
    if not getattr(self.current_program[1].stdout, 'closed', True):
      self.current_program[1].stdout.close()
    if not getattr(self.current_program[1].stderr, 'closed', True): 
      self.current_program[1].stderr.close()
    if len(self.current_program) == 3 and self.current_program[2] is not None:
      from ..misc import LockFile
      LockFile(self.current_program[2].name).remove_stale()
      remove(self.current_program[2].name)
    # now remove reference to current program.
    self.current_program = tuple([None for i in self.current_program])
  def terminate(self):
    """ Terminates current process. """
    if self.current_program[0] is None: return
    try: self.current_program.terminate()
    except: pass
    self._cleanup()

  def kill(self):
    """ Kills current process. """
    if self.current_program[0] is None: return
    try: self.current_program.kill()
    except: pass
    self._cleanup()

  def _poll_current_process(self):
    """ Returns True if no error found.

        If no process is running, there is no error.
    """ 
    if self.current_program[0] is None: return True
    poll = self.current_program[0].poll()
    if poll is None: return
    self._cleanup()
    return poll < 0


class ProgramProcess(Process):
  """ Executes :py:class:`~lada.misc.Program` in child process. """
  def __init__(self, program, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    from copy import deepcopy
    from ..misc import RelativePath
    super(Process, self).__init__(maxtrials, comm, **kwargs)
    self.program = program
    """ Folder to execute. """
    self.current_program = None, None
    """ Tuple specifying currently running process. 

        The tuple contains the following items:
        
        - Popen object.
        - :py:class:`~lada.misc.Program` instance defining the job.
    """ 
    self.poll()

  def poll(): 
    """ Polls current job. """
    from subprocess import Popen
    from shlex import split as split_cmd
    from ..misc import Changedir
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    found_error = self._poll_current_process()
    if self.current_program[0] is not None and not found_error:
      raise Process.StopIteration();

    # increment errors if necessary and check without gone over max trials.
    if found_error: self.nberrors += 1
    if self.nberrors >= self.maxtrials: raise Process.Failed()

    # Open stdout and stderr if necessary.
    with Changedir(program.directory) as cwd:
     file_out = open(program.stdout, "w") if program.stdout is not None else None 
     file_err = open(program.stderr, "w") if program.stderr is not None else None 

    # creates commandline
    d = {}
    d['program'] = self.program.program
    d['cmdline'] = self.program.cmdline 
    d.update(self.comm if self.comm is not None else default_comm)
    d.update(self.params)
    cmdline = split_cmd(mpirun_exe.format(**d))
    process = Popen(cmdline, stdout=file_out, stderr=file_err, cwd=program.directory)
    self.current_program = process, self.program

class UnsharedJobFolderProcess(Process):
  """ Executes folder in child process.
  
      Expects folder with a functional which does not have an iter method.
  """
  def __init__(self, jobfolder, outdir, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    from copy import deepcopy
    from ..misc import RelativePath
    super(UnsharedJobFolderProcess, self).__init__()
    self.current_program = None, None, None
    """ Tuple specifying currently running process. 

        The tuple contains the following items:
        
        - Popen object.
        - :py:class:`~lada.misc.Program` instance defining the job.
        - file to temporary jobfolder current job.
    """ 
    self.comm = comm
    """ MPI communicator. """
    # start first process.
    self.poll()

  def poll(): 
    """ Polls current job. """
    from tempfile import NamedTemporaryFile
    from os.path import join, exists, abspath
    from pickle import dump
    from sys import executable, pypath
    from ..misc import Program
    from . import JobFolder
    from ..lada.jobfolder import __file__ as jobfolder_file
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    found_error = self._poll_current_program()

    # Check directly for errors if possible.
    directory = join(self.outdir, self.jobfolder.name[1:])
    if hasattr(jobfolder.functional, 'Extract') and exists(directory):
      try: found_error = not functional.Extract(directory).success
      except: found_error = True

    # increment errors if necessary and check without gone over max trials.
    if found_error: self.nberrors += 1
    if self.nberrors >= self.maxtrials: raise Process.Failed()

    # create temporary jobfolder file.
    jobfolder = JobFolder() / self.jobfolder.name
    for key, value in self.jobfolder.__dict__.iteritems():
      if key != 'children': jobfolder.__dict__[key] = value
    with NamedTemporaryFile(dir=directory, prefix='.lada_jobtmp', delete=False) as tempfile:
      dump(self.jobfolder)

    # now create process.
    script = abspath(join(basename(jobfolder_file), 'runfolder.py'))
    cmdline = [script, abspath(tempfile.name)] 
    cmdline.extend('--path {0}'.format(path) for path in pypath)
    program = Program(executable, cmdline, directory, None, None)
    process = ProgramProcess(program, comm=self.comm, **self.params)
    self.current_program = process, program, tempfile

class IteratorProcess(Process):
  """ Executes folder in child process.
  
      Expects folder with a functional which has an iter method.
      This method should yield either an extractor object for previously
      successful runs, or a :py:class:`misc.Program` object to execute.
  """
  def __init__(self, jobfolder, outdir, maxtrials=1, comm=None, **kwargs):
    """ Initializes a process. """
    from os.path import join
    from copy import deepcopy
    from ..misc import RelativePath
    super(IteratorProcess, self).__init__(maxtrials, comm, **kwargs)
    self.jobfolder = jobfolder
    """ Folder to execute. """
    self.outdir = RelativePath(join(outdir, self.jobfolder.name[1:])).path
    """ Execution directory of the folder. """
    self.poll()

  def poll(): 
    """ Polls current job. """
    from misc import Program
    # check if we have currently running process.
    # catch StopIteration exception signaling that process finished.
    # does not catch Process.Failed.
    try: self._poll_current_process()
    except Process.StopIteration: pass
    # At this point, loop until find something to do.
    found = False
    params = self.jobfolder.params.copy()
    params.update(self.params)
    for i, program in self.jobfolder.iterator(**params):
      if not getattr(program, 'success', False): 
        found = True
        break;
    # stop if no more jobs.
    if found == False: raise Process.StopIteration()
    # if stopped on or before previous job, then we have a retrial.
    if i <= self.iter_index:
      if self.nberrors >= self.maxtrials: raise Process.Failed()
      self.nberrors += 1

    program = Program(program.program, program.cmdline, program.directory, fileout, filerr)
    process = ProgramProcess(program, comm=self.comm)
    self.current_program = process, program


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
    from misc import Changedir, Program
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
class ExecAll(object):
  """ Executes job-folders on disk. """
  def __init__(self, comm, nbpools=1, retrial=2, folders=None):
    super(ExecAll, self).__init__()
    self.nbpools = nbpools
    """ Number of parallel processes. """
    self.folders = set([])
    """ Path to executable jobfolders. """
    if isintance(folders, list): self.folders = list(folders)
    elif folders is not None: self.folders = [folders]
    self.processes = []
    """ List of currently running processes. """
    self.finished = []
    """ Finished job-folders. """
    self.comm = comm 
    """ Communicator for each job. """
  def add_folder(self, path):
    """ Adds job-folder path to list. """
    from misc import RelativePath
    self.folders.add(RelativePath(path).path)

  def next(self, exhaust=False):
    """ Executes a job-folder if one is still available. """
    from random import choice
    # loop over current jobs.
    for i, (folderpath, jobname, trial, iterator, process) in enumerate(list(self.processes)): 
      poll = process.poll()
      if poll is None: continue
      self.pop(i)
      if poll < 0 and trial < self.trial: # ran into error, retrial.
        self.launch_process(folderpath, jobname, trial+1, iterator)
      else: self.launch_process(folderpath, jobname, trial+1, iterator)

    # add as many jobs as possible.
    if len(folders) == 0: return True
    while len(self.processes) < self.nbpools:
      folderpath = choice(self.folders)
      with findone(folderpath) as job:
        if job != None: self.launch_process(folderpath, job.name, trial, None)
        else: 
          i = self.folders.index(folderpath)
          self.finished.append(self.folderpath.pop(i))
    # If exhaust is True, then returns True if no more processes are running
    # and no more jobs left to execute.
    # If exhaust is False, then returns True if no more jobs left to execute
    # and fewer jobs than pools are left.
    return len(self.processes) == 0 if exhaust else len(self.processes) != self.nbpools
         
  def launch_process(self, folderpath, jobname, trial, iterator):
    """ Adds a process to the list. """
    from ..misc import Program
    if iterator is None:
      with findone(folderpath, jobname) as job:
        if hasattr(job.functional, 'iter'): 
          iterator = job.functional(comm=self.comm, **job.params)
          found = False
          for program in iterator:
            if not getattr(program, 'success', False): 
              found = True
              break
          if not found: return
        else: program = Program(







def execall(path, nbpools=1, processes=None, waittime=2, **kwargs):
  """ Executes job-folders on disk. """
  # loop over executable folders, without risk of other processes detrimentally
  # accessing the file on disk.
  if processes is None: processes = []
  for job in iterall(path):
    if len(processes) >= nbpools:


    
