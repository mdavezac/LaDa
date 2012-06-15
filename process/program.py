from .process import Process
class ProgramProcess(Process):
  """ Executes :py:class:`Program` in child process. """
  def __init__( self, program, outdir, cmdline=None, stdout=None, 
                stderr=None, maxtrials=1, dompi=False, 
                cmdlmodifier=None, **kwargs ):
    """ Initializes a process. """
    from ..error import ValueError
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
    if cmdlmodifier is not None and not hasattr(cmdlmodifier, '__call__'):  
      raise ValueError('cmdlmodifier should None or a callable')
    self.cmdlmodifier = cmdlmodifier
    """ A function to modify command-line parameters.
    
        This function is only invoked for mpirun programs.
        It can be used to, say, make sure a program is launched only with an
        even number of processes. It should add 'placement' to the dictionary.
    """
    self._modcomm = None
    """ An optional modified communicator. 

        Holds communicator optionally returned by commandline communicator.
    """

  def poll(self): 
    """ Polls current job. """
    from . import Fail
    if super(ProgramProcess, self).poll(): return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
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
    if super(ProgramProcess, self).start(comm): return True
    self._next()
    return False
 
  def _next(self):
    """ Starts an actual process. """
    from subprocess import Popen
    from shlex import split as split_cmd
    from ..misc import Changedir
    from ..error import ValueError
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
      if not hasattr(self, '_comm'):
        raise ValueError( "Requested mpi but without passing communicator"     \
                          "(Or communicator was None)." )
      d = {}
      d['program'] = program
      d['cmdline'] = self.cmdline 
      d.update(self._comm)
      if self.cmdlmodifier is not None:
        self._modcomm = self.cmdlmodifier(d, self._comm)
        if self._modcomm is self._comm: self._modcomm = None
      d['placement'] = placement( self._comm if self._modcomm is None          \
                                  else self._modcomm )
      cmdline = split_cmd(mpirun_exe.format(**d))
      cmdline.extend(str(u) for u in self.cmdline)
    else: cmdline = [program] + self.cmdline
    self.process = Popen( cmdline, stdout=file_out,
                          stderr=file_err, cwd=self.outdir )

  def _cleanup(self):
    """ Cleanup files and crap. """
    try: 
      if not getattr(self._stdio[0], 'closed', True):
        self._stdio[0].close()
      if not getattr(self._stdio[1], 'closed', True):
        self._stdio[1].close()
    finally: self._stdio = None, None
    # delete modified communicator, if it exists
    if self._modcomm is not None:
      self._modcomm.cleanup()
      self._modcomm = None
    # general cleanup, including  self._comm
    super(ProgramProcess, self)._cleanup()

  def wait(self):
    """ Waits for process to end, then cleanup. """
    super(ProgramProcess, self).wait()
    self.process.communicate()
    self.poll()

