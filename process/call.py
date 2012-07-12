from .process import Process
class CallProcess(Process):
  """ Calls functional in child python process.
  
      This process pickles_ a callable and its arguments and executes it in a
      child python process.

      .. _pickles: http://docs.python.org/library/pickle.html
  """
  def __init__( self, functional, outdir, stdout=None, stderr=None,\
                maxtrials=1, dompi=False, **kwargs ):
    """ Initializes a process.
    
        :param functional:
          A python callable. It should also be pickle-able.
        :param str outdir: 
          Path where the python child process should be executed.
        :param str stdout:
          Optional path to an output file where the callable's output shall be
          streamed.
        :param str stderr:
          Optional path to an error file where the callable's errors shall be
          streamed.
        :param bool dompi:
          Whether the python child process shouold be launched as an MPI
          process.
        :param int maxtrials:
          Maximum number of times to try re-launching each process upon
          failure. 
        :param kwargs:
          Keyword arguments to the callables should be given here, as keyword
          arguments to :py:class:`CallProcess`.
    """
    from ..misc import RelativePath
    super(CallProcess, self).__init__(maxtrials=maxtrials)
    self.functional = functional
    """ Functional to execute. """
    try: self.functional = self.functional
    except: pass
    self.outdir = RelativePath(outdir)
    """ Execution directory of the folder. """
    self.outdir = RelativePath(outdir).path
    """ Directory where to run job. """
    self.stdout = stdout
    """ Name of standard output file, if any. """
    self.stderr = stderr
    """ Name of standard error file, if any. """
    self.dompi = dompi
    """ Whether to run with mpi or not. """
    self.params = kwargs.copy()
    """ Extra parameters to pass on to iterator. """

  def poll(self):
    """ Polls current job. """
    from . import Fail

    # checks whether program was already started or not.
    if super(CallProcess, self).poll(): return True

    # check whether functional is set. 
    if self.functional is None: return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    try: 
      if self.process.poll() == False: return False
    except Fail:
      self._cleanup()
      if hasattr(self.functional, 'Extract') \
         and self.functional.Extract(self.outdir).success: 
        return True
    else: 
      self._cleanup()
      if not hasattr(self.functional, 'Extract'): return True
      if self.functional.Extract(self.outdir).success: return True
    # at this point, an error must have occured.
    self.nberrors += 1
    if self.nberrors >= self.maxtrials: raise Fail()
    self._next()

  def _next(self):
    """ Launches actual calculation. """
    from pickle import dumps
    from sys import executable, path as pypath
    from tempfile import NamedTemporaryFile
    from .program import ProgramProcess
    from ..misc import Changedir
    # creates temp input script.
    with Changedir(self.outdir) as outdir: pass
    with NamedTemporaryFile(dir=self.outdir, suffix='.py', delete=False) as stdin:
      if self.dompi:
        params = self.params
        stdin.write( "from sys import path\n"\
                     "from boost.mpi import world\n"\
                     "path[:] = {0!r}\n\n"\
                     "from pickle import loads\n"\
                     "params, functional = loads({1!r})\n\n"\
                     "params['comm'] = world\n"\
                     "functional(**params)\n"\
                     .format(pypath, dumps((params, self.functional))) )
      else: 
        params = {'comm': self._comm}
        params.update(self.params)
        stdin.write( "from sys import path\n"\
                     "path[:] = {0!r}\n\n"\
                     "from pickle import loads\n"\
                     "params, functional = loads({1!r})\n\n"\
                     "functional(**params)\n"\
                     .format(pypath, dumps((params, self.functional))) )
      self._stdin = stdin.name

    # now create process. maxtrials is one if Extract exists, so that we can
    # check success using that instead.
    self.process = ProgramProcess( executable, cmdline=[self._stdin], 
                                   outdir=self.outdir, stdout=self.stdout,
                                   stderr=self.stderr, maxtrials=1,
                                   dompi=self.dompi )
    self.process.start(comm=self._comm.lend('all') if self.dompi else None)
    return False

  def start(self, comm):
    """ Starts current job. """
    if super(CallProcess, self).start(comm): return True
    # Check directly for errors/success if possible.
    if hasattr(self.functional, 'Extract') \
       and self.functional.Extract(self.outdir).success \
       and self.process is None \
       and self.params.get('overwrite', False) == False:
      return True
    self._next()
    return False;

  def _cleanup(self):
    """ Removes temporary script. """
    from os import remove
    super(CallProcess, self)._cleanup()
    if not hasattr(self, '_stdin'): return
    try: remove(self._stdin)
    except: pass
    finally: del self._stdin

  def wait(self):
    """ Waits for process to end, then cleanup. """
    if super(CallProcess, self).wait(): return True
    while self.poll() == False: self.process.wait()
    self._cleanup()
    return False

