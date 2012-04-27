from .process import Process
class CallProcess(Process):
  """ Calls functional in child python process. """
  def __init__( self, functional, outdir, cmdline=None, stdout=None, stderr=None,\
                maxtrials=1, dompi=False, **kwargs ):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(CallProcess, self).__init__(maxtrials=maxtrials, **kwargs)
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

  def poll(self, wait=False): 
    """ Polls current job. """
    from . import Fail

    # checks whether program was already started or not.
    if super(CallProcess, self).poll(): return True

    # check whether functional is set. 
    if self.functional is None: return True
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    try: 
      if wait == True: self.process.wait()
      elif self.process.poll() == False: return False
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
                                   nompi=self.dompi )
    self.process.start(comm=self._comm if self.dompi else None)
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
    while not self.poll(): self.process.wait()
    return False

