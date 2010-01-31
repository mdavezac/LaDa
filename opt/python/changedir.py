""" Class to change working directory within local context only.
   
    >>> with Changedir(path) as pwd:
    >>>   ...
"""
class Changedir:
  """ Works with "with" statement to temporarily change the working directory 
   
      >>> with Changedir(path) as pwd:
      >>>   ...
  """
  def __init__(self, pwd): self.pwd = pwd

  def __enter__(self):
    """ Changes working directory """
    from os import getcwd, chdir
    from os.path import exists, isdir
    
    self.oldpwd = getcwd()

    assert exists(self.pwd) and isdir(self.pwd), "Could not find working directory."
    chdir(self.pwd)

    return self.pwd

  def __exit__(self, type, value, traceback):
    """ Moves back to old pwd """
    from os import chdir
    from os.path import exists, isdir

    assert exists(self.oldpwd) and isdir(self.oldpwd), "Old directory does not exist anymore."
    chdir(self.oldpwd)
