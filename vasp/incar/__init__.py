""" Subpackage defining vasp incar parameters. """
__docformat__ = "restructuredtext en"
__all__ = [ "SpecialVaspParam", "NElect", "Algo", "Precision", "Ediff",\
            "Encut", "FFTGrid", "Restart", "UParams", "IniWave", 'Ediffg', "EncutGW", \
            "Incar", "Magmom", 'Npar', 'Boolean', 'Integer', 'Choices', 'PrecFock' ]
from _params import SpecialVaspParam, NElect, Algo, Precision, Ediff,\
                    Encut, FFTGrid, Restart, UParams, IniWave, Magmom,\
                    Npar, Boolean, Integer, PrecFock, NonScf, Ediffg, Choices, \
                    EncutGW
from ...misc import add_setter


class Incar(object):
  """ Contains vasp Incar parameters. 

      The following assumes you know how to write an INCAR. Although you won't
      need to anymore.  This class separates vasp parameters from methods to
      launch and control vasp.

      There are two kinds of parameters: 

        - Normal parameters which will simply print "NAME = VALUE" to the incar
        - Special parameters which enhance the default behavior of vasp
      
      The special parameters are always called first. They may change the
      values of one or more "normal" vasp parameters. For instance
      :py:attr:`restart <Incar.restart>` will set :py:attr:`istart
      <Incar.istart>` and :py:attr:`icharg <Incar.icharg>` depending on
      the availability of the wavefunctions and charge densitie.
  """
  def __init__(self): 
    super(Incar, self).__init__()

    # first, actually sets these two variables by hand, since they are used in __setattr__.
    super(Incar, self).__setattr__("params", {})
    super(Incar, self).__setattr__("special", {})
    self.add_param = "ispin",       1 
    self.add_param = "isif",        1
    self.add_param = "ismear",      None
    self.add_param = "isigma",      None
    self.add_param = "nsw",         None
    self.add_param = "ibrion",      None
    self.add_param = "potim",       None
    self.add_param = "nbands",      None
    self.add_param = "lorbit",      None
    self.add_param = "addgrid",     None
    self.add_param = "isym",        None
    self.add_param = "symprec",     None
    self.add_param = "nupdown",     None
    self.add_param = "lmaxmix",     4
    self.add_param = "lmaxfockae",  None
    self.add_param = "nomega",      None
    self.add_param = "istart",      None
    self.add_param = "icharge",     None
    # objects derived from SpecialVaspParams will be recognized as such and can
    # be added without further fuss.
    self.nelect      = NElect(0)
    self.algo        = Algo()
    self.precision   = Precision("accurate")
    self.ediff       = Ediff(1e-4)
    self.ediffg      = Ediffg(None)
    self.encut       = Encut(None)
    self.encutgw     = EncutGW(None)
    self.fftgrid     = FFTGrid(None)
    self.restart     = Restart(None)
    self.U_verbosity = UParams("occupancy")
    self.magmom      = Magmom()
    self.npar        = Npar(None)
    self.precfock    = PrecFock(None)
    self.nonscf      = NonScf(False)

    self.lwave       = Boolean("lwave", False)
    self.lcharg      = Boolean("lcharg", True)
    self.lvtot       = Boolean("lvtot", False)
    self.lrpa        = Boolean("lrpa", None)
    self.loptics     = Boolean("loptics", None)
    self.lpead       = Boolean("lpead", None)
    self.nelm        = Integer("nelm", None)
    self.nelmin      = Integer("nelmin", None)
    self.nelmdl      = Integer("nelmdl", None)

  def incar_lines(self, **kwargs):
    """ List of incar lines. """

    structure = kwargs['structure']
    # prints system name. This is not an option!
    result = []
    if hasattr(structure, "name"):
      if len(structure.name) != 0:
        result.append("SYSTEM = \"{0}\"\n".format(structure.name))

    # gathers special parameters.
    # Calls them first in case they change normal key/value pairs.
    specials = []
    for key, value in self.special.items():
      if value.value is None: continue
      line = value.incar_string(**kwargs)
      if line is not None: specials.append(line + "\n")
    # prints key/value pairs
    for key, value in self.params.items():
      if value is None: continue
      if isinstance(value, bool):  value = ".TRUE." if value else ".FALSE."
      else: 
        try: value = str(value)
        except ValueError: 
          raise ValueError("Could not convert vasp parameter {0} to string: {1}.".format(key, value))
      result.append( "{0: >18s} = {1}\n".format(key.upper(), value))
    # adds special parameter lines.
    result.extend(specials)
    return result

  @add_setter
  def add_param(self, args):
    """ Adds/sets a vasp parameter.
    
        Consists of a key value pair. 

        >>> vasp.add_param = "ispin", 2

        This will result in the INCAR as "ISPIN = 2". Once set, the value can be accessed directly:

        >>> vasp.add_param = "ispin", 2
        >>> vasp.ispin = 1
        >>> print vasp.ispin # prints 1
    """
    key, value = args
    if isinstance(value, SpecialVaspParam):
      if key in self.params: del self.params[key] # one or other dictionary.
      self.special[key] = value
    else:
      if key in self.special: del self.special[key] # one or other dictionary.
      self.params[key] = value

  def __getattr__(self, name): 
    """ Gets a VASP parameter from standard and special dictionaries. """
    if name in self.params: return self.params[name]
    elif name in self.special: return self.special[name].value
    raise AttributeError("Unknown parameter " + name)

  def __setattr__(self, name, value):
    """ Sets a VASP parameter to standard and special dictionaries. """
    if isinstance(value, SpecialVaspParam):
      if name in self.params: del self.params[name]
      self.special[name] = value
    elif name in self.params: self.params[name] = value
    elif name in self.special: self.special[name].value = value
    else: super(Incar, self).__setattr__(name, value)

  def __delattr__(self, name): 
    """ Deletes a VASP parameter from standard and special dictionaries. """
    if name in self.__dict__: return self.__dict__.pop(name)
    elif name in self.params: return self.params.pop(name)
    elif name in self.params: return self.special.pop(name).value
    raise AttributeError("Unknown vasp attribute " + name + ".")

  def __dir__(self):
    result = [u for u in self.__dict__ if u[0] != '_'] 
    result.extend([u for u in self.params.keys() if u[0] != '_'])
    result.extend([u for u in self.special.keys() if u[0] != '_'])
    return list(set(result))

  @add_setter
  def symmetries(self, value):
    """ Type of symmetry used in the calculation.
  
        This sets :py:attr:`isym` and :py:attr:`symprec` vasp tags. Can be
        "off" or a float corresponding to the tolerance used to determine
        symmetry operation. 
    """
    if value is None: self.isym = None
    elif str(value).lower() == "off" or value is "0" or value is False: self.params["isym"] = 0
    elif str(value).lower() == "on" or value is True or value is True:
       self.symprec = None
       self.isym = None
    elif isinstance(value, float): 
       self.symprec = value
       self.isym = None
    else: raise ValueError("Uknown value when setting symmetries ({0}).".format(value))

  @add_setter
  def smearing(self, args):
    """ Value of the smearing used in the calculation. 
  
        It can be specified as a string:
          
        >>> vasp.smearing = "type", x
       
        Where type is any of "fermi", "gaussian", "mp", "tetra", "metal", or "insulator",
        and x is the energy scale in eV.

        - fermi use a Fermi-Dirac broadening.
        - gaussian uses a gaussian smearing.
        - mp is for Methfessel-Paxton, it should be specified as "mp N x",
          where N is the order the mp method.
        - tetra tetrahedron method without Bloechl correction.
        - bloechl means tetrahedron method with Bloechl correction.
        - "metal x" is equivalent to "mp 1 x"
        - insulator is equivalent to "tetra bloechl".
        - if x is omitted a default value of 0.2eV is used.
    """
    if args is None: 
      self.ismear, self.isigma = None, None
      return

    if isinstance(args, str): args = args.split()
    first = args[0]
    has_second = len(args) > 1
    second = args[1] if len(args) > 1 else None
    has_third = len(args) > 2
    third = args[2] if len(args) > 2 else None

    if first is None:
      self.ismear = None
      if has_second: self.isigma = None
    elif first == "fermi" or first == "-1":    
      self.ismear = -1
      if has_second: self.isigma = second
    elif first == "gaussian" or first == "0":
      self.ismear = 0
      if has_second: self.isigma = second
    elif first == "metal":
      self.ismear = 1
      if has_second: self.isigma = second
    elif first == "mp" or first == "metal":
      if has_third:
        self.ismear = second
        self.isigma = third
      elif has_second: self.isigma = second
      else: self.ismear = 1
      assert self.ismear >= 1, "Mehtfessel-Paxton order must be at least 1."
    elif first == "bloechl" or first == "-5" or first == "insulator":
      self.ismear = -5
      if has_second: self.isigma = second
    elif first == "tetra" or first == "-4":
      self.ismear = -4
      if has_second: self.isigma = second
    else: 
      try: self._value = int(first)
      except: raise ValueError("Unknown smearing value {0}.\n".format(value))
      if self._value < 1:
        raise ValueError("Unknown smearing value {0}.\n".format(value))

  @property
  def relaxation(self):
    """ Sets type of relaxation.
    
        It can be set to a single value, or to a tuple of up to four elements:

        >>> vasp.relaxation = "static" 
        >>> vasp.relaxation = "static", 20
      
        - first argument can be "static", or a combination of "ionic",
          "cellshape", and "volume". The combination must be allowed by
          `ISIF`__.
        - second (optional) argument is `nsw`_
        - third (optional) argument is `ibrion`_
        - fourth (optional) argument is `potim`_

        .. note:: Some combinations will raise an error, eg asking for
           ionic relaxation but setting nsw to zero. However, if you are
           feeling creative and know what you are doing, you can always set the
           parameters by hand:

           >>> vasp.nsw, vasp.isif = 0, 1

        .. __: http://cms.mpi.univie.ac.at/vasp/guide/node112.html
        .. _nsw: http://cms.mpi.univie.ac.at/vasp/guide/node108.html
        .. _ibrion: http://cms.mpi.univie.ac.at/vasp/guide/node110.html
        .. _potim: http://cms.mpi.univie.ac.at/vasp/vasp/POTIM_tag.html
    """
    nsw = 0 if self.nsw is None else self.nsw
    if self.ibrion is None: ibrion = -1 if nsw < 0 else 0
    else: ibrion = self.ibrion
    if self.isif is None: isif = 0 if ibrion == 0 else 2
    else: isif = self.isif
   
    if nsw < 2 or ibrion == -1: return "static"
    result = ""
    if isif < 5: result += "ionic "
    if isif > 2 and isif < 7: result += "cellshape "
    if isif in [3, 6, 7]: result += "volume"
    return result.lstrip().rstrip()

  @relaxation.setter
  def relaxation(self, args): 
    """ Sets the kind of relaxation to be performed. """
    import re

    dof =  args.lower() if isinstance(args,str) else str(args[0]).lower()
    ionic = re.search( "ion(ic|s)?", dof ) is not None
    cellshape = re.search( "cell(\s+|-|_)?(?:shape)?", dof ) is not None
    volume = re.search( "volume", dof ) is not None

    nsw, ibrion, potim = None, None, None
    if not isinstance(args, str):
      if len(args) > 1: nsw = int(args[1])
      if len(args) > 2: ibrion = int(args[2])
      if len(args) > 3: potim = int(args[3])

    # static calculation.
    if (not ionic) and (not cellshape) and (not volume):
      if dof != 'static':
        raise RuntimeError("Unkown value for relaxation: {0}.".format(arg))
      self.params["isif"] = 1
      self.params["ibrion"] = -1
      if ibrion is not None and ibrion != -1:
        raise ValueError("Cannot set ibrion to anything but -1 for static calculations.")
      if nsw is not None and nsw != 0:
        raise ValueError("static calculation with nsw > 0 is way too creative.")
      self.params["nsw"] = None
      if potim is not None: self.params["potim"] = potim

    else: # Some kind of relaxations. 
      # ionic calculation.
      if ionic and (not cellshape) and (not volume):   self.params["isif"] = 1
      elif ionic and cellshape and (not volume):       self.params["isif"] = 4
      elif ionic and cellshape and volume:             self.params["isif"] = 3
      elif (not ionic) and cellshape and volume:       self.params["isif"] = 6
      elif (not ionic) and cellshape and (not volume): self.params["isif"] = 5
      elif (not ionic) and (not cellshape) and volume: self.params["isif"] = 7
      elif ionic and (not cellshape) and volume: 
        raise RuntimeError, "VASP does not allow relaxation of atomic position"\
                            "and volume at constant cell-shape.\n"

      if ibrion is None and self.params["ibrion"] in [None, -1]: self.params["ibrion"] = 2
      elif ibrion is not None: 
        if ibrion == -1: raise ValueError("Cannot set ibrion to -1 with strain relaxations.")
        if ibrion == 0 and self.params["isif"] == 1: 
          raise ValueError("Cannot set ibrion to 0 with strain relaxations.")
        self.params["ibrion"] = ibrion
      if nsw is not None:
        if nsw <= 0: raise ValueError("Cannot set nsw < 1 and perform strain relaxations.")
        self.params["nsw"] = nsw
      elif self.params["nsw"] is None or self.params["nsw"] == 0: self.params["nsw"] = 50
      if potim is not None: self.params["potim"] = potim
      if self.ediffg is not None:
        if self.ediffg < self.ediff and self.ediffg > 0: self.ediffg = None

  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("params")
    special = d.pop("special")
    return d, params, special
  def __setstate__(self, args):
    super(Incar, self).__setattr__("params", args[1])
    super(Incar, self).__setattr__("special", args[2])
    d = self.__dict__.update(args[0])
