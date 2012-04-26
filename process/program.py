from .process import Process
class ProgramProcess(Process):
  """ Executes :py:class:`Program` in child process. """
  def __init__( self, program, outdir, cmdline=None, stdout=None, 
                stderr=None, maxtrials=1, dompi=False, **kwargs ):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(ProgramProcess, self).__init__(maxtrials, **kwargs)
    self.program = program
    """ External program to execute. """
    self.outdir = RelativePath(outdir).path
    """ Directory where to run job. """
    self.cmdline = [] if cmdline is None else cmdline
    """ Command line for the program. """
    self.stdout = stdout
    """ Name of standard output file, if any. """
    self.stderr = stderr
    """ Name of standard error file, if any. """
    self.dompi = dompi
    """ Whether to run with mpi or not. """
    self._stdio = None, None
    """ Standard output/error files. """

  def poll(self): 
    """ Polls current job. """
    from . import Fail
    if super(ProgramProcess, self).poll(): return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    if self.process is not None:
      poll = self.process.poll()
      if poll is None: return False
      elif poll != 0:
        self.nberrors += 1
        if self.nberrors >= self.maxtrials:
          self._cleanup()
          raise Fail(poll)
      else:
        self._cleanup()
        return True

    # increment errors if necessary and check without gone over max trials.
    self._next() # restart process.
    return False

  def start(self, comm):
    """ Starts current job. """
    if super(ProgramProcess, self).start(comm): return 
 
  def _next(self):
    """ Starts an actual process. """
    from subprocess import Popen
    from shlex import split as split_cmd
    from ..misc import Changedir
    from .. import mpirun_exe, placement
    from . import which
    # Open stdout and stderr if necessary.
    with Changedir(self.outdir) as cwd:
      file_out = None if self.stdout is None else open(self.stdout, "w") 
      file_err = None if self.stderr is None else open(self.stderr, "w") 
      self._stdio = file_out, file_err

    # creates commandline
    program = which(self.program)
    if self.dompi: 
      d = {}
      d['program'] = program
      d['cmdline'] = self.cmdline 
      d.update(self._comm)
      d.update(self.params)
      d['placement'] = placement(self._comm)
      print mpirun_exe.format(**d)
      cmdline = split_cmd(mpirun_exe.format(**d))
      cmdline.extend(str(u) for u in self.cmdline)
    else: cmdline = [program] + self.cmdline
    self.process = Popen(cmdline, stdout=file_out, stderr=file_err, cwd=self.outdir)

  def _cleanup(self):
    """ Cleanup files and crap. """
    try: 
      if not getattr(self._stdio[0], 'closed', True):
        self._stdio[0].close()
      if not getattr(self._stdio[1], 'closed', True):
        self._stdio[1].close()
    finally: self._stdio = None, None
    super(ProgramProcess, self)._cleanup()

  def wait(self):
    """ Waits for process to end, then cleanup. """
    self.process.communicate()
    self.poll()

