""" Subpackage defining vasp incar parameters. """
from _params import Standard, NoPrintStandard, AlgoValue, \
                    PrecValue, EdiffValue, EncutValue, \
                    SmearingValue, SymValue, FFTValue, \
                    RestartValue, RelaxationValue

class Incar(object):
  """ Contains vasp Incar parameters. 

      This class separates vasp parameters from methods to launch and control
      vasp. vasp attributes can be listed by iterating over this class, or
      calling iter.
  """ 

  def __init__(self): 
    object.__init__(self) # calls base class.

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
    self.precision = PrecValue()
    """ Sets accuracy of calculation. 
        Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
    """
    self.prec = EdiffValue()
    """ Sets the convergence criteria for electronic minimization.

        This tolerance is divided by the number of atoms in the system.  For
        this reason, printing to incar is doned via the return to __call__.
    """
    self.nsw = Standard("NSW", 0, validity = lambda x: int(x) == float(x) and int(x) >= 0)
    """ Maximum number of ionic steps. \n """ 
    self.encut = EncutValue(safety=1.25)
    """ Gets maximum cutoff from POTCAR.
        Actual printed value is that times the safety. 
    """
    self.smearing = SmearingValue()
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
    self.fftgrid = FFTValue(grid = None)
    """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """
    self.relaxation = RelaxationValue()
    """ Sets ISIF in incar depending on type relaxation required. 
    
          - if set to None or empty string, then no relaxation.
          - if ionic is in string, then includes ionic relaxation.
          - if cellshape is in string, then includes cell relaxation.
          - if volume is in string, then includes volume relaxation.
    
        Makes sure that the parameters make sense together. 
        Can also be set using an integer between 0 and 7. See VASP manual. 
    """

  restart = RestartValue(None)
  """
      Directory where to restart, or None.
      
      If None, then starts from scratch.
      If this directory contains WAVECAR or CHGCAR, then restarts from
      wavefunctions or charge. If this directory does not exist, issue an
      error.
  """

  def __iter__(self):
    """ Iterates over vasp incar parameters.

        To be identified as an incar parameter, an attribute should have an
        incar_string function. In this way, the Incar class can be subclassed
        while retaining this method.
    """

    for name in dir(self):
      if name[0] == '_': continue
      attr = getattr(self, name)
      if hasattr(attr, "incar_string"): yield attr
