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
  def __init__(self, workdir = None, comm = None, keep = False):
    self.workdir = workdir
    self.keep = keep
    self.comm = comm

  def __enter__(self):
    """ Creates temporary directory """
    from os import mkdir, listdir, makedirs
    from os.path import exists, isdir
    from tempfile import mkdtemp
    from shutil import rmtree 
    from boost.mpi import broadcast

    is_root = self.comm == None
    if self.comm != None: is_root = self.comm.rank == 0
    if self.workdir != None:
      if not exists(self.workdir): makedirs(self.workdir)
      assert exists(self.workdir) and isdir(self.workdir),\
             "Could not create working directory."
    if is_root: self._tempdir = mkdtemp(dir=self.workdir)
    else: self._tempdir  = None
    if self.comm != None:
      assert broadcast(self.comm, "Tempdir: am in sync", root = 0) == "Tempdir: am in sync", \
             "Processes not in sync"
      self._tempdir = broadcast(self.comm, value=self._tempdir, root=0)
      for i in range(1, self.comm.size):
        self.comm.barrier()
        if i != self.comm.rank: continue
        if not exists(self._tempdir): mkdir(self._tempdir)
    assert exists(self._tempdir) and isdir(self._tempdir),\
           "Could not create temporary working directory."
    assert len(listdir(self._tempdir)) == 0,\
           "Could not create temporary working directory."
    return self._tempdir

  def __exit__(self, type, value, traceback):
    """ Deletes temporary directory """
    from shutil import rmtree
    from os.path import exists, isdir
    if self.keep: return
    if self.comm != None:
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
