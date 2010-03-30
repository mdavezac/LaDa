""" Defines specie specific methods and objects. """


def U(type = 1, l=2, U=0e0, J=0e0 ):
  """ Creates an LDA+U parameter 

      LDA+U is always LSDA+U here. 
      @param type: (1|2|"liechtenstein"|"dudarev")
      @param l: (0|1|2|"s"|"p"|"d") channel for which to apply U. Default: 2.
  """

  if hasattr(type, "lower"):
    type = type.lower()
    if type == "liechtenstein": type = 1
    elif type == "dudarev": type == 2
  if hasattr(l, "lower"):
    l = l.lower()
    if l == "s": l = 0
    elif l == "p": l == 1
    elif l == "d": l == 2
  try: l = int(l)
  except: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  if l < 1 or l > 2: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  return { "type": int(type), "l": l, "U": U, "J": J, "func": "U" }

def nlep(type = 1, l=2, U=0e0, Eref=None ):
  """ Creates nlep parameters 

      LDA+U is always LSDA+U here. 
      @param type: (1|2|"liechtenstein"|"dudarev")
      @param Eref: (float|None) reference energy for E-nlep, if specified.
        Default:None.
      @param l: (0|1|2|"s"|"p"|"d") channel for which to apply U. Default: 2.
  """
  if hasattr(type, "lower"):
    if type.lower() == "liechtenstein": type = 1
    elif type.lower() == "dudarev": type == 2
  if Eref == None: 
    return { "type": int(type), "l": l, "U": U, "func": "nlep" }
  if hasattr(l, "lower"):
    l = l.lower()
    if l == "s": l = 0
    elif l == "p": l == 1
    elif l == "d": l == 2
  try: l = int(l)
  except: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  if l < 1 or l > 2: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  else: 
    return { "type": int(type), "l": l, "U": U, "ref": Eref, "func": "enlep" }


class Specie(object):
  """ Holds atomic specie information:  """
  def __init__(self, symbol, path, U=None):
    """ Initializes a specie.
        @param symbol: is the atomic symbol
        @type symbol: str
        @param path: to the directory with the potcar for this particular atomic types.
          This directory should contain a POTCAR or POTCAR.Z file.
        @type path: str
        @param U: LDA+U parameters. It should a list of dictionaries, one entry
                  per momentum channel. The simplest approach is to use
                  >>> Species(..., U=[U(l=1, U=1.5), U(l=2, U=2.5)],...)
    """
    import os.path

    self.symbol = symbol
    self.path = os.path.expanduser( path )
    if U == None: self.U = []
    else: self.U = U



  @property
  def enmax(self):
    """ Maximum recommended cutoff """
    import re
    from os.path import exists
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
    from os.path import exists
    if not exists(os.path.join(self.path, "POTCAR")):
      raise IOError, "Could not find potcar in " + self.path
    with open(os.path.join(self.path, "POTCAR"), "r") as potcar:
      file.readline()
      return float(file.readline().split()[0]) # shoud be number on second line
