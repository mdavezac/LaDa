""" Class to handle creation/destruction of temporary directories. 
   
    >>> with Tempdir(path) as tempdir:
    >>>   ...
"""
class Tempdir:
  """ Works with "with" statement to create/destroy a temporary directory 
   
      >>> with Tempdir(path) as tempdir:
      >>>   ...
      @param workdir: working directory where to create a directory. 
      @param keep: when True, does not destroy directory on exit, eg for debug
         purposes.
  """
  def __init__(self, workdir = None, comm = None, keep = False, debug = None):
    from ..mpi import Communicator
    self.workdir = workdir
    self.keep = keep
    self.comm = Communicator(comm)
    self.debug = debug

  def __enter__(self):
    """ Creates temporary directory """
    from os import mkdir, listdir, makedirs
    from os.path import exists, isdir, join
    from tempfile import mkdtemp
    from shutil import rmtree 

    if self.comm.is_root:
      if self.workdir != None:
        if not exists(self.workdir): makedirs(self.workdir)
        assert exists(self.workdir) and isdir(self.workdir),\
               "Could not create working directory."
      if self.debug == None: 
        self._tempdir = mkdtemp(dir=self.workdir)
      else:
        self._tempdir = join(self.workdir, self.debug)
        if not exists(self._tempdir): makedirs(self._tempdir)
    else: self._tempdir  = None
    if self.comm.is_mpi:
      self._tempdir = self.comm.broadcast(self._tempdir)
      for i in range(1, self.comm.size):
        self.comm.barrier()
        if i != self.comm.rank: continue
        if not exists(self._tempdir): mkdir(self._tempdir)
    assert exists(self._tempdir) and isdir(self._tempdir),\
           "Could not create temporary working directory."
    assert self.debug != None or len(listdir(self._tempdir)) == 0,\
           "Could not create temporary working directory."
    self.comm.barrier()
    return self._tempdir

  def __exit__(self, type, value, traceback):
    """ Deletes temporary directory """
    from shutil import rmtree
    from os.path import exists, isdir
    if self.keep: return
    if self.comm.is_mpi:
      for i in range(self.comm.size):
        self.comm.barrier()
        if not i == self.comm.rank: continue
        if not exists(self._tempdir): continue
        if not isdir(self._tempdir): continue
        rmtree(self._tempdir)
    else:
      if not exists(self._tempdir): return
      if not isdir(self._tempdir): return
      rmtree(self._tempdir)
