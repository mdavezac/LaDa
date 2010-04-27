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
    self.pwd = pwd
    self.comm = comm

  def __enter__(self):
    """ Changes working directory """
    from os import getcwd, chdir
    from os.path import exists, isdir
    
    self.oldpwd = getcwd()

    is_root = self.comm == None
    if not is_root: is_root = self.comm.rank == 0
    if is_root: # creates directory if does not exist.
      if not exists(self.pwd): makedirs(self.pwd)
    assert exists(self.pwd) and isdir(self.pwd), "Could not find working directory."
    chdir(self.pwd)

    return self.pwd

  def __exit__(self, type, value, traceback):
    """ Moves back to old pwd """
    from os import chdir
    from os.path import exists, isdir

    assert exists(self.oldpwd) and isdir(self.oldpwd), "Old directory does not exist anymore."
    chdir(self.oldpwd)
