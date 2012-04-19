class CallProcess(Process):
  """ Calls functional in child python process. """
  def __init__( self, functional, outdir, cmdline=None, stdout=None, stdin=None\
                maxtrials=1, comm=None, dompi=False, **kwargs ):
    """ Initializes a process. """
    from ..misc import RelativePath
    super(ProgramProcess, self).__init__(maxtrials=maxtrials, comm=comm, **kwargs)
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
    sef._tempfile = None


  def poll(self): 
    """ Polls current job. """
    from . import Program
    from sys import executable, pypath
    from pickle import dumps
    from .process import Stop, Fail
    from .program import ProgramProcess
    # check whether functional is set. 
    if self.functional is None: raise Stop()
    # check if we have currently running process.
    # if current process is finished running, closes stdout and stdout.
    found_error = not self._poll_process()
    # Check directly for errors if possible.
    if hasattr(jobfolder.functional, 'Extract') and exists(self.outdir):
      try: found_error = not functional.Extract(self.outdir).success
      except: found_error = True

    # increment errors if necessary and check without gone over max trials.
    if found_error: self.nberrors += 1
    if self.nberrors >= self.maxtrials: raise Fail()

    if self.dompi: params = self.params
    else: 
      params = {'comm': self.comm}
      params.update(self.params)
    with NamedTemporaryFile(dir=self.outdir, prefix='.lada_script', delete=False) as tempfile:
      tempfile.write( "from sys import path\n"\
                      "path[:] = {0!r}\n\n"\
                      "from pickle import loads\n"\
                      "params, functional = loads({1!r})\n\n"\
                      "functional(*params)\n".format(path, dumps((params, self.functional))) )
      self._tempfile.name = tempfile.name

    # now create process.
    self.process = ProgramProcess(executable, self.outdir, cmdline=[sef._tempfile], 
                                  stdout=self.stdout, stderr=self.stderr, maxtrials=1,
                                  comm=self.comm if self.dompi else None, 
                                  nompi=self.dompi)
    self.process.start()

  def start(self):
    """ Starts current job. """
    self.poll()

  def _cleanup(self):
    """ Cleans up after process. """
    from os import remove
    if self._tempfile is not None:
      try: remove(self._tempfile)
      except: pass
      finally: self._tempfile = None
    super(ProgramProcess, self)._cleanup()
