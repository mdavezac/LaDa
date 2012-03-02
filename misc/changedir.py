""" Class to change working directory within local context only.
   
    >>> with Changedir(path) as pwd:
    >>>   ...
"""
class Changedir:
  """ Works with "with" statement to temporarily change the working directory 
   
      >>> with Changedir(path) as pwd:
      >>>   ...
  """
  def __init__(self, pwd): 
    self.pwd = pwd

  def __enter__(self):
    """ Changes working directory """
    from os import getcwd, chdir, makedirs
    from os.path import exists, isdir
    
    self.oldpwd = getcwd()

    if not exists(self.pwd): makedirs(self.pwd)
    if not exists(self.pwd):
      raise IOError("Could not create working directory {0}".format(self.pwd))
    if not isdir(self.pwd):
      raise IOError("{0} is not a directory.".format(self.pwd))
    chdir(self.pwd)

    return self.pwd

  def __exit__(self, type, value, traceback):
    """ Moves back to old pwd """
    from os import chdir
    from os.path import exists, isdir
    if not (exists(self.oldpwd) or not isdir(self.oldpwd)):
      raise IOError("Old directory does not exist anymore.")
    chdir(self.oldpwd)
