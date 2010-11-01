""" Defines specie specific methods and objects. """
__docformat__ = "restructuredtext en"


def U(type = 1, l=2, U=0e0, J=0e0 ):
  """ Creates an LDA+U parameter 

      :Parameters:
        type : 1|2|"liechtenstein"|"dudarev"
          A string or integer specifying the type of the Hubbard U. Defaults
          to 1.
        l : 0|1|2|"s"|"p"|"d"
          Channel for which to apply U. Defaults to 2.
        U : float
          Hubbard U. Defaults to 0.
        J : float
          Hubbard J. Defaults to 0.

      LDA+U is always LSDA+U here. 
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
  """ Creates `nlep`_ parameters 

      :Parameters:
        type : 1|2|"liechtenstein"|"dudarev"
          A string or integer specifying the type of the Hubbard U. Defaults
          to 1.
        l : 0|1|2|"s"|"p"|"d"
          Channel for which to apply U. Defaults to 2.
        U0 : float
          First nlep parameter. Defaults to 0.
        U1 : float
          Second (e)nlep parameter. Defaults to 0.

      `Non Local External Potentials`__ attempt to correct in part for the
      band-gap problem. It was proposed in PRB **77**, 241201(R) (2008).
      :note: LDA+U is always LSDA+U here. 

      .. _nlep : http://dx.doi.org/10.1103/PhysRevB.77.241201
      __ nlep_
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
  """ Holds atomic specie information.
  
      Instances of this object define an atomic specie for VASP calculations.
      In addition, it may contain information used to build a set of
      high-throughput jobs.
  """
  def __init__(self, directory, U=None, oxidation=None, **kwargs):
    """ Initializes a specie.

        :Parameters:
          directory : str
            Directory with the potcar for this particular atomic types.  This
            directory should contain an *unzipped* POTCAR file.
          U 
            LDA+U parameters. It should a list of dictionaries, one entry per
            momentum channel, and each entry returned by a call to `U`
            or `nlep`.
          oxidation : int
            Maximum oxidation state (or minimum if negative).
          kwargs : dict
            Any other keyworkd argument is added as an attribute of this object.
    """
    from ..opt import RelativeDirectory

    self._directory = RelativeDirectory(directory)
    if oxidation != None: self.oxidation = oxidation
    if U == None: self.U = []
    elif isinstance(U, dict): self.U = [U]
    else: self.U = [u for u in U] # takes care of any kind of iterator.

    # sets up other arguments.
    for k, v in kwargs.items(): setattr(self, k, v)

  @property
  def directory(self):
    """ Directory where the POTCAR file may be found. """
    return self._directory.path
  @directory.setter
  def directory(self, value): self._directory.path = value

  @property 
  def path(self):
    """ Path to POTCAR file. """
    from os.path import join
    return join(self.directory, "POTCAR")

  @property
  def enmax(self):
    """ Maximum recommended cutoff """
    import re
    self.potcar_exists()
    with self.read_potcar() as potcar:
      r = re.compile("ENMAX\s+=\s+(\S+);\s+ENMIN")
      p = r.search(potcar.read())
      if p == None: raise AssertionError, "Could not retrieve ENMAX from " + self.directory
      return float( p.group(1) )
  
  @property
  def valence(self):
    """ Number of valence electrons specified by pseudo-potential """ 
    import re
    self.potcar_exists()
    with self.read_potcar() as potcar:
      potcar.readline()
      return float(potcar.readline().split()[0]) # shoud be number on second line

  def potcar_exists(self):
    """ Raises IOError if POTCAR file does not exist. """
    from os.path import exists
    assert exists(self.path), IOError("Could not find POTCAR in {0}.".format(self.directory))

  def read_potcar(self):
    """ Returns handle/context to POTCAR file. """
    self.potcar_exists()
    return open(self.path, "r") 

  def __repr__(self):
    """ Represents a specie. """
    string = "{0.__class__.__name__}('{0._directory.unexpanded}'".format(self)
    for k, v in self.__dict__.items():
      if k[0] == '_directory': continue
      if k == 'U' and len(v) == 0: continue
      try: assert repr(v)[0] != '<' 
      except: continue
      string += ", {name}={value}".format(name=k, value=repr(v))
    return string + ')'



