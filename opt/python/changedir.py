""" Class to change working directory within local context only.
   
    >>> with Changedir(path) as pwd:
    >>>   ...
"""
class Changedir:
  """ Works with "with" statement to temporarily change the working directory 
   
      >>> with Changedir(path) as pwd:
      >>>   ...
  """
  def __init__(self, pwd, comm = None): 
    from ..mpi import Communicator
    self.pwd = pwd
    self.comm = Communicator(comm)

  def __enter__(self):
    """ Changes working directory """
    from os import getcwd, chdir, makedirs
    from os.path import exists, isdir
    
    self.oldpwd = getcwd()

    if self.comm.is_root and not exists(self.pwd): makedirs(self.pwd)
    self.comm.barrier()
    assert exists(self.pwd), "Could not find working directory %s." % (self.pwd)
    assert isdir(self.pwd), "%s is not a directory." % (self.pwd)
    chdir(self.pwd)
    self.comm.barrier()

    return self.pwd

  def __exit__(self, type, value, traceback):
    """ Moves back to old pwd """
    from os import chdir
    from os.path import exists, isdir

    assert exists(self.oldpwd) and isdir(self.oldpwd), "Old directory does not exist anymore."
    chdir(self.oldpwd)
    self.comm.barrier()
