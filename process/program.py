from .process import Process
class ProgramProcess(Process):
  """ Executes :py:class:`Program` in child process. """
  def __init__( self, program, outdir, cmdline=None, stdout=None, 
                stderr=None, stdin=None, maxtrials=1, dompi=False, 
                cmdlmodifier=None, onfinish=None, **kwargs ):
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
    self.stdin = stdin
    """ Name of standard error file, if any. """
    self.dompi = dompi
    """ Whether to run with mpi or not. """
    self._stdio = None, None, None
    """ Standard output/error/input files. """
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
    self.onfinish = onfinish
    """ Callback when the processes finishes. 

        Called even on error. Should take two arguments:
        
          - process: holds this instance
          - error: True if an error occured.

        It is called before the :py:meth:`_cleanup` method. In other words, the
        process is passed as it is when the error is found.
    """ 

  def poll(self): 
    """ Polls current job. """
    from . import Fail
    if super(ProgramProcess, self).poll(): return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    poll = self.process.poll()
    if poll is None: return False
    # call callback.
    if hasattr(self.onfinish, '__call__'):
      self.onfinish(process=self, error=(poll!=0))
    # now check possible error.
    if poll != 0:
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

  def start(self, comm=None):
    """ Starts current job. """
    if super(ProgramProcess, self).start(comm): return True
    self._next()
    return False
 
  def _next(self):
    """ Starts an actual process. """
    from os import environ
    from ..misc import Changedir
    from ..error import ValueError
    from .. import mpirun_exe, launch_program as launch
    from . import which
    # Open stdout and stderr if necessary.
    with Changedir(self.outdir) as cwd:
      file_out = None if self.stdout is None else open(self.stdout, "w") 
      file_err = None if self.stderr is None else open(self.stderr, "w") 
      file_in  = None if self.stdin is None else open(self.stdin, "r") 
      self._stdio = file_out, file_err, file_in

    # creates commandline
    program = which(self.program)
    if self.dompi: 
      if not hasattr(self, '_comm'):
        raise ValueError( "Requested mpi but without passing communicator"     \
                          "(Or communicator was None)." )
      formatter = {}
      cmdl = ' '.join(str(u) for u in self.cmdline)
      formatter['program'] = '{0} {1}'.format(program, cmdl)
      # gives opportunity to modify the communicator before launching a
      # particular program.
      if self.cmdlmodifier is not None:
        self._modcomm = self.cmdlmodifier(formatter, self._comm)
        if self._modcomm is self._comm: self._modcomm = None
      comm    = self._comm if self._modcomm is None else self._modcomm 
      cmdline = mpirun_exe
    else:
      cmdline   = program + ' '.join(str(u) for u in self.cmdline)
      comm      = None
      formatter = None

    self.process = launch( cmdline, comm=comm, formatter=formatter,
                           env=environ, stdout=file_out, stderr=file_err,
			   stdin=file_in, outdir=self.outdir )

  def _cleanup(self):
    """ Cleanup files and crap. """
    try: 
      if not getattr(self._stdio[0], 'closed', True):
        self._stdio[0].close()
      if not getattr(self._stdio[1], 'closed', True):
        self._stdio[1].close()
      if not getattr(self._stdio[2], 'closed', True):
        self._stdio[2].close()
    finally: self._stdio = None, None, None
    # delete modified communicator, if it exists
    if self._modcomm is not None:
      self._modcomm.cleanup()
      self._modcomm = None
    # general cleanup, including  self._comm
    super(ProgramProcess, self)._cleanup()

  def wait(self):
    """ Waits for process to end, then cleanup. """
    super(ProgramProcess, self).wait()
    self.process.wait()
    self.poll()

