from .process import Process
class ProgramProcess(Process):
  """ Executes :py:class:`Program` in child process. """
  def __init__( self, program, outdir, cmdline=None, stdout=None, stdin=None\
                maxtrials=1, comm=None, dompi=False, **kwargs ):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(Process, self).__init__(maxtrials, comm, **kwargs)
    self.program = program
    """ External program to execute. """
    self.outdir = RelativePath(outdir).path
    """ Directory where to run job. """
    self.cmdline = [] if None else cmdline
    """ Command line for the program. """
    self.stdout = stdout
    """ Name of standard output file, if any. """
    self.stderr = stderr
    """ Name of standard error file, if any. """
    slef.dompi = dompi
    """ Whether to run with mpi or not. """
    self._stdio = None, None
    """ Standard output/input files. """

  def poll(self): 
    """ Polls current job. """
    from subprocess import Popen
    from shlex import split as split_cmd
    from ..misc import Changedir
    from . import Stop, Fail
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    found_error = self._poll_process()
    if self.process is not None and not found_error: raise Stop()

    # increment errors if necessary and check without gone over max trials.
    if found_error: self.nberrors += 1
    if self.nberrors >= self.maxtrials: raise Fail()
    self.start() # restart process.

  def start(self):
    """ Starts current job. """
    # Open stdout and stderr if necessary.
    with Changedir(self.outdir) as cwd:
      file_out = open(self.program.stdout, "w")\
                    if self.program.stdout is not None else None 
      file_err = open(self.program.stderr, "w")\
                    if self.program.stderr is not None else None 
      self._stdio = file_out, file_err

    # creates commandline
    if self.dompi: cmdline = self.cmdline
      d = {}
      d['program'] = self.program.program
      d['cmdline'] = self.program.cmdline 
      d.update(self.comm if self.comm is not None else default_comm)
      d.update(self.params)
      cmdline = split_cmd(mpirun_exe.format(**d))
    else: cmdline = self.cmdline
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
