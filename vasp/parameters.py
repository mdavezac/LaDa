from parameter_types import Standard, NoPrintStandard, AlgoValue, PrecValue, EdiffValue

class Incar(object):
  """ Contains vasp Incar parameters. 
      This class separates vasp parameters from methods to launch and control
      vasp. vasp attributes can be listed by iterating over this class, or
      calling iter.
  """ 
  _exclude_from_iteration = []
  """ Public attributes over which not to iterate. 
      This means these public attributes are not proper vasp parameters. """

  def __iter__(self):

    vasp = [ u for u in dir(self) if u[0] != '_' ] # removes private stuff.
    vasp = [ u for u in vasp if u not in self._exclude_from_iteration ]
    def generator():
      for i in vasp: yield (getattr(self, i), i)
    return generator()
        

 
  def __init__(self):

    self.iniwave = Standard("INIWAV", "random", validity = lambda x: x=="random" or x=="jellium")
    """ Initializes wave functions with \"random\"(default) or \"jellium\" """ 
    self.nelect = NoPrintStandard("NELECT", 0 , validity = lambda x: x >= 0)
    """ Sets number of electrons in calculation.
        0(default) lets VASP compute it from species in the system. 
    """
    self.nbands = NoPrintStandard("NBANDS", 0, validity = lambda x: x >= 0)
    """ Sets number of bands to include in calculation.
        0(default) lets VASP decide.
    """ 
    self.potim = Standard( "POTIM", 0.5, validity = lambda x: float(x) > 0e0)
    """ Sets ionic motion time step. """ 
    self.nspins = Standard( "ISPIN", 1, 
                            validity = lambda x: int(x) == float(x) and (int(x)==1 or int(x)==2) )
    """ Sets number of spins. Must be either 1 or 2. """
    self.algo = AlgoValue()
    """ Electronic minimization. 
        Can be \"very fast\", \"fast\", or \"normal\" (default). 
    """ 
    self.value = PrecValue()
    """ Sets accuracy of calculation. 
        Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
    """
    self.prec = EdiffValue(self)
    """ Sets the convergence criteria for electronic minimization.
        This tolerance is divided by the number of atoms in the system. 
        For this reason, printing to incar is doned via the return to __call__.\n"
    """
    self.nsw = Standard("NSW", 0, validity = lambda x: int(x) == float(x) and int(x) >= 0)
    """ Maximum number of ionic steps. \n """ 
    self.encut = EcutValue(self, safety=1.25)
    """ Gets maximum cutoff from POTCAR.
        Actual printed value is that times the safety. 
    """
    self.smearing = SmearingValue(string = "insulator" )
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
    self.isym = SymValue()
    """ Type of symmetry used in the calculation.
        Can be "off" or a float corresponding to the tolerance used to determine
        symmetry operation.  By default, it is 1e-5.
    """
    self.fft = FFTValue(self, grid = None)
    """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """

