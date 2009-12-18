class Incar(object):
  """ Contains vasp Incar parameters. 
      This class separates vasp parameters from methods to launch and control
      vasp. Any attribute here should be printable as a vasp incar.
  """

  class Standard(object):
    """ A standard parameter in the form of key/value pair. """
    def __init__(self, key, value, validity=None, doc=None):

      self.key = key
      self.value = value
      if doc != None: self.__doc__ = doc

    def __set__(self, parent, value):
      """ Sets value of key/value pair """
      if self.validity != None:
        assert self.validity(value), \
               "Value %s = %s is invalid.\n%s\n"  % (self.key, value, self.__doc__)
      self.value = value

    def __get__(self, parent, parenttype=None):
      """ Gets value of key/value pair """
      return self.value

    def __str__(self):
      """ Prints in INCAR format. """
      return "%s = %s" % (self.key, self.value)

    def add_to_doc(self, parent):

      lock = "_LockDocOfObject%s" % (self.__name__)
      if hasattr(parent.__class__, lock): return 
      setattr(parent, lock, None)
      parent.__doc__ += "self.%s:\n%s" % (self.__name__, self.%s)

  def NoPrintStandard(Standard):
    """ Does not print if value is the default given on initialization. """
    def __init__(self, key, value, validity=None, doc=None):
      Standard.__init__(self, key, value, validity, doc)
      self.default = value
    def __str__(self):
      if self.default == self.value: 
        return "# %s = VASP default." % (self.key)
      return Standard.__str__(self)



  def __init__(self):

    self.iniwave = Standard( "INIWAV", "random", validity = lambda x: x=="random" or x=="jellium",\
                             doc = """Initializes wave functions with \"random\" or \"jellium\"""" ) 
#   self.nelect = NoPrintStandard( "NELECT", 0, validity = lambda x: x >= 0,
#                                  """ Sets number of electrons in calculation.
#                                      0 lets VASP compute it from species in the system. 
#                                  """ )
#   self.nbands = NoPrintStandard( "NBANDS", 0, validity = lambda x: x >= 0,
#                                  """ Sets number of bands to include in calculation.
#                                      0 lets VASP decide.
#                                  """ )
#   self.potim = Standard( "POTIM", 0.5, validity = lambda x: float(x) > 0e0,
#                          """ Sets ionic motion time step. """ )
#   self.nspins = Standard( "ISPIN", 1, validity = lambda x: int(x)==1 or int(x)=2,\
#                           """ Sets number of spins. Must be either 1 or 2. """ )

#   # defines a class for possible electronic minimization algorithms
#   class AlgoValue(Standard): 
#     """ Electronic minimization. 
#         Can be \"very fast\", \"fast\", or \"normal\" (default). 
#     """ 
#     def __init__(self): Standard.__init__( self, "ALGO", "Normal") 
#     def __set__(self, parent, value):
#       lower = value.lower()
#       lower = lower.replace("_", " ")
#       if lower == "very fast": self.value = "Very_Fast"
#       elif lower == "fast": self.value = "Fast"
#       elif lower == "normal": self.value = "Normal"
#       else: raise ValueError, "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)
#   self.algo = AlgoValue()

#   # defines a class for possible algorithms
#   class PrecValue(Standard):
#     """ Sets accuracy of calculation. 
#         Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
#     """
#     def __init__(self): Standard.__init__( self, "PREC", "accurate") 
#     def __set__(self, parent, value):
#       self.value = value.lower()
#       assert    self.value == "accurate" \
#              or self.value == "low" \
#              or self.value == "medium" \
#              or self.value == "high",
#              "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)
#   self.value = PrecValue()

#   # Defines the tolerance of the electronic minimization. 
#   class EdiffValue(Standard):
#     """ Sets the convergence criteria for electronic minimization.
#         This tolerance is divided by the number of atoms in the system. 
#     """
#     def __init__(self, parent):
#       Standard.__init__("PREC", 1e-4, validity = lambda x: x > 0e0)
#       self.parent = parent
#     def __str__(self):
#       return "PREC = %f " % (self.value / float(len(self.parent.system.atoms)))
#   self.prec = EdiffValue(self)



#   self.encut = -1
#   self.cutoff_safety = 1.25
#   self.nelect  = 0
#   self.nbands = 0
#   self.nspins = 1
#   self.smearing = "gaussian 0.2"
#   self.kpoint = "d=10 o=0.5"
#   self.potim = 0.5
#   self.isym  = "default 1e-5"

#   self.algo = "fast"
#   self.prec = "accurate"
#   self.ediff = 1e-4
#   self.fft = (-1,-1,-1)

#   self.relaxation = "global volume ionic cellshape"
#   self.nsw = 40

#   # directory where infiles are stored.
#   self.indir = ""

#   self.other = {}
