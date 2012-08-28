from .process import Process
class ProgramProcess(Process):
  """ Executes an external program

      This process creates and manages the execution of an external program,
      say VASP_ or CRYSTAL_, via a `subprocess.Popen`__ instance. The external
      program can be launched with or without MPI, with or without standard
      output/input/error files. It is always launched within a specified
      directory. 
  
      A typical use case, taken from process/test/program.py, is the
      following:

      .. code-block :: python

        program = ProgramProcess( executable, outdir=dir, 
                                  cmdline=['--sleep', 0, '--order', 4], 
                                  stdout=stdout, dompi=True )
        program.start(comm)
        try: program.wait()
        except Fail:
          # do something

      The above launches an external program taking a set of arguments. Its
      output is piped  to a specific file for later grepping. It is launched
      using the super-computer's MPI interface, with the number of processors
      specified by ``comm``. It is launched in a directory ``dir``. 
      The snippet above has python wait for the external program to finish,
      while checking for exceptions if the program fails.
      The external program is started only once :py:meth:`start` is called. 

      .. note ::
        
        A :py:class:`~lada.process.Fail` exception is thrown when the program
        returns with a non-zero exit code. However, some propietary MPI
        crapware, such as Cray's. will return 0 whenever ``MPI::Finalize()`` is
        called, even when the program itself returns non-zero. As a result, it
        is not possible to rely on a :py:class:`~lada.process.Fail` exception
        being thrown correctly on all machines at all times.

      .. __ : http://docs.python.org/library/subprocess.html#subprocess.Popen
      .. _VASP: http://www.vasp.at/
      .. _CRYSTAL: http://www.crystal.unito.it/
  """
  def __init__( self, program, outdir, cmdline=None, stdout=None, 
                stderr=None, stdin=None, maxtrials=1, dompi=False, 
                cmdlmodifier=None, onfinish=None, onfail=None, **kwargs ):
    """ Initializes a process.
    
        :param str program: 
          Path to the executable of interest.
        :param str outdir:
          Path to the directory where the program should be executed.
        :param list cmdline: 
          List of commandline arguments. The elements of the list should be
          translatable to strings in a meaningfull way, via str_.
        :param str stdout:
          Path to a file where the standard output of the executable should go.
          Ignored if None, in which the standard output is likely piped to
          `sys.stdout <http://docs.python.org/library/sys.html#sys.stdout>`_.
        :param str stderr: 
          Path to a file where the standard error of the executable should go.
          Ignored if None, in which the standard error is likely piped to
          `sys.stderr <http://docs.python.org/library/sys.html#sys.stderr>`_.
        :param str stdin: 
          Path to a file where the standard input of the executable should go.
          Ignored if None, in which the standard input is likely piped to
          `sys.stdin <http://docs.python.org/library/sys.html#sys.stdind>`_.
          Note that there is no way to communicate with the process managed by
          this instance *via* the standard input. Indeed,
          `subprocess.Popen.communicate
          <http://docs.python.org/library/subprocess.html#subprocess.Popen.communicate>`_
          locks the process until it finishes, which defeat the whole purpose
          of the :py:mod:`process` module.
        :param int maxtrials: 
          Maximum number of restarts. If the external program fails, it will
          automically be restarted up to this number of times. 
        :param bool dompi:
          If True, then the external program is launched with MPI.
          The MPI infrastructure should be set up correctly, meaning
          :py:func:`lada.launch_program`,
          :py:func:`lada.machine_dependent_call_modifier`,
          :py:func:`lada.modify_global_comm`, :py:data:`lada.mpirun_exe`,
          :py:data:`~lada.default_comm`, :py:data:`~lada.figure_out_machines`.
        :param cmdlmodifier:
          This function is called prior to launching the program. It can be
          used to modify the formatting dictionary. It should return a
          communicator to use in the actual call. This communicator can be the
          same as passed to the function. However, if the communicator needs to
          be modified, then a new one should be created (and the original left
          untouched). Ownership should not be retaine beyond the call to this
          function.
        :type cmdlmodifier: cmdlmodifier(formatter, comm)->comm
        :param onfinish:
          Called once the process is finished. The first argument is this
          instance, the second is True if an error occurred. 
        :type onfinish: onfinish(self, exitcode!=0)->None
        :param kwargs:
          Other keyword arguments are ignored.

        .. _str : http://docs.python.org/library/functions.html#str
    """
    from os.path import getsize
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
    self.onfail = onfail
    """ Called if program fails. 

	Some program, such as CRYSTAL, return error codes when unconverged.
	However, does not necessarily mean the program failed to run. This
	function is called when a failure occurs, to make sure it is real or
	not. It should raise Fail if an error has occurred and return normally
        otherwise.
    """ 

  def poll(self): 
    """ Polls current job.
    
        :returns: True if external program is finished.
        :raises Fail: If external program returns a non-zero exit code.
    """
    from . import Fail
    if super(ProgramProcess, self).poll(): return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    poll = self.process.poll()
    if poll is None: return False
    # call callback.
    if self.onfinish is not None:
      try: self.onfinish(process=self, error=(poll!=0))
      except Exception as e: 
        import sys, traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        tb = traceback.format_tb(exc_traceback)
        self._cleanup()
        raise Fail( 'Error on call to "onfinish"\n{0}: {1}\n{2}\n'             \
                    .format(type(e), e, '\n'.join(tb)) )
    # now check possible error.
    if poll != 0:
      if self.onfail is not None:
        try: self.onfail(process=self, error=poll)
        except Fail: 
          if self.nberrors >= self.maxtrials:
            self._cleanup()
            raise
          else: self.nberrors += 1
        else:
          self._cleanup()
          return True
      else:
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
    if super(ProgramProcess, self).start(comm): return True
    self._next()
    return False
  start.__doc__ = Process.start.__doc__
 
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
      cmdl = ' '.join(str(u) for u in self.cmdline)
      cmdline   = '{0} {1}'.format(program, cmdl)
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

