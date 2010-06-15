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
    elif type == "dudarev": type = 2
  if hasattr(l, "lower"):
    l = l.lower()
    assert len(l) == 1, "Uknown input %s." % (l)
    if   l[0] == 's': l = 0
    elif l[0] == 'p': l = 1
    elif l[0] == 'd': l = 2
  try: l = int(l)
  except: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  if l < 1 or l > 2: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  return { "type": int(type), "l": l, "U": U, "J": J, "func": "U" }

def nlep(type = 1, l=2, U0=0e0, U1=None ):
  """ Creates nlep parameters 

      LDA+U is always LSDA+U here. 
      @param type: (1|2|"liechtenstein"|"dudarev")
      @param U0: (float) first E-nlep parameter. Defaults to 0e0.
      @param U1: (float|None) second E-nlep parameter, if specified. Otherwise
                 reverts to nlep. Default:None.
      @param l: (0|1|2|"s"|"p"|"d") channel for which to apply U. Default: 2.
  """
  if hasattr(type, "lower"):
    type = type.lower()
    if type == "liechtenstein": type = 1
    elif type == "dudarev": type = 2
  if hasattr(l, "lower"):
    l = l.lower()
    assert len(l) == 1, "Uknown input %s." % (l)
    if   l[0] == 's': l = 0
    elif l[0] == 'p': l = 1
    elif l[0] == 'd': l = 2
  try: l = int(l)
  except: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  if l < 0 or l > 2: raise ValueError, "Moment l should be 0|1|2|s|p|d." 
  elif U1 == None: 
    return { "type": int(type), "l": l, "U": U0, "func": "nlep" }
  else: 
    return { "type": int(type), "l": l, "U0": U0, "U1": U1, "func": "enlep" }


class Specie(object):
  """ Holds atomic specie information:  """
  def __init__(self, symbol, path, U=None, oxidation=None):
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
    self.path = os.path.abspath(os.path.expanduser(path))
    if oxidation != None: self.oxidation = oxidation
    if U == None: self.U = []
    elif isinstance(U, dict): self.U = [U]
    else: self.U = U



  @property
  def enmax(self):
    """ Maximum recommended cutoff """
    import re
    from os.path import exists, join
    if not exists(join(self.path, "POTCAR")):
      raise IOError, "Could not find potcar in " + self.path
    with open(join(self.path, "POTCAR"), "r") as potcar:
      r = re.compile("ENMAX\s+=\s+(\S+);\s+ENMIN")
      p = r.search(potcar.read())
      if p == None: raise AssertionError, "Could not retrieve ENMAX from " + self.path
      return float( p.group(1) )
  
  @property
  def valence(self):
    """ Number of valence electrons specified by pseudo-potential """ 
    import re
    from os.path import exists, join
    if not exists(join(self.path, "POTCAR")):
      raise IOError, "Could not find potcar in " + self.path
    with open(join(self.path, "POTCAR"), "r") as potcar:
      potcar.readline()
      return float(potcar.readline().split()[0]) # shoud be number on second line

  def __str__(self):
    """ Prints atomic symbol only. """
    return str(self.symbol)
