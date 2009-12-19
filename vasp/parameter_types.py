""" Paremeters types for use as attributes in Incar """
class Standard(object):
  """ A standard parameter in the form of key/value pair. """
  def __init__(self, key, value, validity=None):

    self.key = key
    self.value = value

  def __set__(self, owner, value):
    """ Sets value of key/value pair """
    if self.validity != None:
      assert self.validity(value), \
             "Value %s = %s is invalid.\n%s\n"  % (self.key, value, self.__doc__)
    self.value = value

  def __get__(self, owner, parenttype=None):
    """ Gets value of key/value pair """
    return self.value

  def __str__(self):
    """ Prints in INCAR format. """
    return "%s = %s" % (self.key, self.value)

class NoPrintStandard(Standard):
  """ Does not print if value is the default given on initialization. """
  def __init__(self, key, value, validity=None):
    Standard.__init__(self, key, value, validity)
    self.default = value
  def __str__(self):
    if self.default == self.value: 
      return "# %s = VASP default." % (self.key)
    return Standard.__str__(self)

class AlgoValue(Standard): 
  """ Electronic minimization. 
      Can be \"very fast\", \"fast\", or \"normal\" (default). 
  """ 
  def __init__(self): Standard.__init__( self, "ALGO", "Normal") 
  def __set__(self, owner, value):
    lower = value.lower()
    lower = lower.replace("_", " ")
    if lower == "very fast": self.value = "Very_Fast"
    elif lower == "fast": self.value = "Fast"
    elif lower == "normal": self.value = "Normal"
    else: raise ValueError, "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)

class PrecValue(Standard):
  """ Sets accuracy of calculation. 
      Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
  """
  def __init__(self): Standard.__init__( self, "PREC", "accurate") 
  def __set__(self, owneowner, value):
    self.value = value.lower()
    assert    self.value == "accurate" \
           or self.value == "low"      \
           or self.value == "medium"   \
           or self.value == "high",    \
           "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)

class EdiffValue(Standard):
  """ Sets the convergence criteria for electronic minimization.
      This tolerance is divided by the number of atoms in the system. 
      For this reason, printing to incar is doned via the return to __call__.\n"
  """
  def __init__(self, owner):
    Standard.__init__(self, "PREC", 1e-4, validity = lambda x: x > 0e0)
    self._owner = owner
  def __str__(self):
    return "PREC = %f " % (self.value / float(len(self._owner.system.atoms)))


class EncutValue(Object):
  """ Defines energy cutoff for calculation. 
      The base cutoff is obtained from the species of the owner specified in
      __init__. The actual cutoff is that times the safety (also specified in __init__).
  """
  class Safety(object):
    def __init__(self, value):
      self.value = value
    def __set__(self, owner, value):
      assert value > 0, "Ecut safety must be larger than 0.\n"
      self._value = value
    def __get__(self, owner, parenttype=None): return self._value

  def __init__(self, owner, safety = 1.25):
    self.safety = Safety(safety)
    self._owner = owner

  @property
  def value(self):
    """ Retrieves max ENMAX from list of species. """
    assert hasattr(self, "species"), "No species member in owner."
    import os.path
    import subprocess
    import re
    from math import ceil

    result = 0
    for s in self._owner.species: 
      stdout = ""
      if os.path.exists( os.path.join(s.path, "POTCAR") ):
        filename = os.path.join(s.path, "POTCAR")
        cmd = subprocess.Popen(["grep", "ENMAX", s.path + "POTCAR"], \
                               stdout=subprocess.PIPE)
        stdout = cmd.stdout.readline()
      elif os.path.exists( os.path.join(s.path, "POTCAR.Z") ):
        filename = os.path.join(s.path, "POTCAR.Z")
        cmd0 = subprocess.Popen(["zcat", filename], \
                                stdout=subprocess.PIPE)
        cmd = subprocess.Popen(["grep", "ENMAX"], \
                               stdin=cmd0.stdout, stdout=subprocess.PIPE)
        stdout = cmd.communicate()[0]
      else: raise AssertionError, "Could not find potcar in " + s.path
  
      r = re.compile("ENMAX\s+=\s+(\S+);\s+ENMIN")
      p = r.search(stdout)
      if p == None: raise AssertionError, "Could not retrieve ENMAX from " + s.path
      if result < float( p.group(1) ): result = float( p.group(1) )

    return ceil(result)

  def __str__(self):
    """ Prints as INCAR input. """
    print "ENCUT = %f" % (self.value*self.safety)


class SmearingValue(object):
  """ Value of the smearing used in the calculation. 
      It can be specified as a string: "type x", where type is any of fermi,
      gaussian, mp, tetra, metal or insulator, and x is the energy scale in eV.
          - fermi use a Fermi-Dirac broadening.
          - gaussian uses a gaussian smearing.
          - mp is for Methfessel-Paxton, it should be specified as "mp N x",
            where N is the order the mp method.
          - tetra tetrahedron method without Bloechl correction.
          - bloechl means tetrahedron method with Bloechl correction.
          - metal is equivalent to "mp 1"
          - insulator is equivalent to "tetra bloechl".
          - if x is omitted a default value of 0.2eV is used.
  """
  def __init__(self, string = "insulator"):
    self.type = _Type(string)
    self.value = _Value(string)

  def __set__(self, owner, value):

    first = None
    second = None 
    third = None

    if isinstance(value, str):
      data = value.split()
      first = data[0].lower()
      if len(data) > 1: second = data[1]
      if len(data) > 2: third  = data[2]
      if len(data) > 3: raise "Unknown input %s. smearing: %s\n" % (value, self.__doc__)
    elif hasattr(value, __len__):
      first = str(value[0]).lower()
      if len(value) > 1: second = value[1]
      if len(value) > 2: third  = value[2]
      if len(value) > 3: raise "Unknown input %s. smearing: %s\n" % (value, self.__doc__)
    elif isinstance(value, float):
      self.value = float(value)
      return
    elif isinstance(value, int): raise "int is not an accepted type for smearing."

    if first == "fermi" or first == "-1":    
      self.type == -1
      if second != None: self.value = float(second)
    elif first == "gaussian" or first == "0":
      self.type == 0
      if second != None: self.value = float(second)
    elif first == "mp":
      if third != None:
        self.type = int(second)
        self.value = float(third)
      elif second != None:
        self.type = int(second)
      else: self.type = 1
      assert self.type >= 1, "Mehtfessel-Paxton order must be at least 1."
    elif first == "bloechl" or first == "-5":
      self.type = -5
      if second != None: self.value = float(second)
    elif first == "tetra" or first == "-4":
      self.type = -4
      if second != None: self.value = float(second)
    else: 
      try: self._value = int(first)
      except: raise "Unknown smearing value %s.\nsmearing: %s\n" % (value, self.__doc__)
      assert self._value >= 1, "Unknown smearing value %s.\nsmearing: %s\n" % (value, self.__doc__)

  def __str__(self):
    return "ISMEAR  = %s\nSIGMA   = %f\n" % (self.type, self.value)

  
class SymValue(object):
  """ Type of symmetry used in the calculation.
      Can be "off" or a float corresponding to the tolerance used to determine
      symmetry operation.  By default, it is 1e-5.
  """
  def __init__(self, string = 1e-5)
    self.__set__(string)

  def __set__(self, owner, value):

    strval = str(value)
    if strval == "off": self._isym = 0
    elif:
      self._isym = None
      self._tol = float(value) 
  def __str__(self):
    if self._isym == None:
      return "SYMPREC = %3.2e" % (self._tol)
    return "ISYM = 0\nSYMPREC = %3.2e" % (self._tol)
  
class FFFTValue(object):
    """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """
    def __init__(self, owner, grid = None):
      self._grid = grid
      self._owner = _owner 

    def __set__(self, owner, value):
      self._owner = owner
      self._grid = value
      if self._grid == None: return
      assert len(self._grid) == 3: "FFT grid should be none or a 3-tuple.\n"
      assert self._grid[0] > 0: "FFT grid components should be positive.\n"
      assert self._grid[1] > 0: "FFT grid components should be positive.\n"
      assert self._grid[2] > 0: "FFT grid components should be positive.\n"

    def __str__(self):
      grid = self._grid 
      if grid == None: raise "Not implemented."
      return "NGX = %i\n NGY = %i\n NGZ = %i" % tuple(grid)


  

