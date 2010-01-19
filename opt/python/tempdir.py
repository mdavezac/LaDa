""" Class to handle creation/destruction of temporary directories. 
   
    >>> with Tempdir(path) as tempdir:
    >>>   ...
"""
class Tempdir:
  """ Works with "with" statement to create/destroy a temporary directory 
   
      >>> with Tempdir(path) as tempdir:
      >>>   ...
  """
  def __init__(self, workdir): self.workdir = workdir

  def __enter__(self):
    """ Creates temporary directory """
    from os import makedirs
    from os.path import exists, isdir
    from tempfile import mkdtemp
    from shutil import rmtree 

    if not exists(self.workdir): makedirs(self.workdir)
    assert exists(self.workdir) and isdir(self.workdir),\
           "Could not create working directory."
    self._tempdir = mkdtemp(dir=self.workdir)
    assert exists(self._tempdir) and isdir(self._tempdir),\
           "Could not create temporary working directory."
    return self._tempdir

  def __exit__(self, type, value, traceback):
    """ Deletes temporary directory """
    from shutil import rmtree
    if exists(self._tempdir) and isdir(self._tempdir):
      rmtree(self._tempfile)
