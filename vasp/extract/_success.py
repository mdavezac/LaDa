""" Checks for sucess of vasp calculation """
class Success(object):
  """ Checks for success of vasp calculation """
  def __init__(self): object.__init__(self)

  def __get__(self, owner, type=None):
    from os.path import exists, join
    from .. import files
    import re

    for path in [files.OUTCAR, files.CONTCAR]:
      if owner.directory != "": path = join(owner.directory, path)
      if not exists(path): return False
      
    path = files.OUTCAR 
    if len(owner.directory): path = join(owner.directory, path)

    with open(path, "r") as file:
      regex = re.compile(r"""General\s+timing\s+and\s+accounting
                             \s+informations\s+for\s+this\s+job""", re.X)
      for line in file:
        if regex.search(line) != None: return True
    return False

  def __set__(self, owner, value):
    raise TypeError, "Success cannot be set. It is a result only.\n"
