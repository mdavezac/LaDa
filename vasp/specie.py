""" Defines specie methods. """
class Specie:
  """ Holds atomic specie information:  """
  def __init__(self, symbol, path ):
    """ Initializes a specie.
        @param symbol: is the atomic symbol
        @type symbol: str
        @param path: to the directory with the potcar for this particular atomic types.
          This directory should contain a POTCAR or POTCAR.Z file.
        @type path: str
    """
    import os.path

    self.symbol = symbol
    self.path = os.path.expanduser( path )

  @property
  def enmax(self):
    """ Maximum recommended cutoff """
    import re
    import os.path import exists
    if not exists(os.path.join(self.path, "POTCAR")):
      raise IOError, "Could not find potcar in " + self.path
    with open(os.path.join(self.path, "POTCAR"), "r") as potcar:
      r = re.compile("ENMAX\s+=\s+(\S+);\s+ENMIN")
      p = r.search(potcar.read())
      if p == None: raise AssertionError, "Could not retrieve ENMAX from " + self.path
      return float( p.group(1) )

  @property
  def valence(self):
    """ Number of valence electrons specified by pseudo-potential """ 
    import re
    import os.path import exists
    if not exists(os.path.join(self.path, "POTCAR")):
      raise IOError, "Could not find potcar in " + self.path
    with open(os.path.join(self.path, "POTCAR"), "r") as potcar:
      file.readline()
      return float(file.readline().split()[0]) # shoud be number on second line
