""" Subpackage defining vasp incar parameters. """
__docformat__ = "restructuredtext en"
__all__ = [ "SpecialVaspParam", "NElect", "Algo", "Precision", "Ediff",\
            "Encut", "FFTGrid", "Restart", "UParams", "IniWave",\
            "Incar", "Magmom", 'Npar', 'Boolean', 'Integer', 'Choices',\
            'PartialRestart', 'PrecFock', 'NonScf' ]
from _params import SpecialVaspParam, NElect, Algo, Precision, Ediff,\
                    Encut, FFTGrid, Restart, UParams, IniWave, Magmom,\
                    Npar, Boolean, Integer, PrecFock, NonScf, PartialRestart
from ...opt.decorators import add_setter


class Incar(object):
  """ Contains vasp Incar parameters. 

      The following assumes you know how to write an INCAR. Although you won't
      need to anymore.  This class separates vasp parameters from methods to
      launch and control vasp. vasp attributes can be listed by iterating over
      this class, or calling iter.

      Parameters(attributes) which are present by default are the following:
         - ``ispin``: Sets number of spins. Must be either 1 or 2. 
         - ``ismear``: Smearing function. Can be set with property `set_smearing`. 
         - ``sigma``: Smearing parameter. Can be set with property `set_smearing`.
         - ``isif``: Degrees of freedom to relax. Can be set using `relaxation`. 
         - ``nsw``: Number of ionic steps. Can be set using `relaxation`. 
         - ``ibrion``: ionic-relaxation method. Can be set using `relaxation`. 
         - ``potim``: ionic-relaxation step. Can be set using `relaxation`. 
         - ``iniwave``: initial wavefunction to use can be either "random" or "jellium".   
         - ``nelect``: sets number of electrons in calculation above and beyond valence.
             - 0(default) lets VASP compute it from species in the system. 
             - 1 would mean +1 electron
             - -1 would mean +1 hole
             - etc
         - ``algo``: electronic minimization. Can be \"very fast\", \"fast\" (default), or \"normal\".
         - ``precision``: sets accuracy of calculation. Can be \"accurate\"
           (default), \"low\", \"medium\", \"high\".
         - ``ediff``: sets the tolerance per atom of electronic minimization.
            This tolerance is multiplied by the number of atoms in the system,
            eg consistent from one system to another.
         - ``encut``: factor by which ENMAX of species is multiplied.
         - ``fftgrid``: Either a 3-tuple of integers setting NGX and friend, or
           anything else (other than None). In the latter case, a fake VASP
           calculation is performed to get VASP recommended values.
         - ``symprec``: tolerance when determining symmetries.
         - ``magmom``: Sets magnetic moment. See `incar.Magmom`.
         - ``npar``: Sets npar.
         - ``lcorr``: Defaults to None.
         - ``lplane``: Defaults to None.
         - ``nbands``: Defaults to None.
         - ``lorbit``: Defaults to None.
         - ``lsorbit``: Defaults to None.
         - ``addgrid``: Defaults to None.
         - ``nupdown``: Defaults to None.
         - ``loptics``: Defaults to None.
         - ``lmaxmix``: Defaults to 4.
         - ``magmom``: Defaults to None.
         - ``weimin``: Defaults to 0.
         - ``lwave`` : Defaults to None. 
         - ``lmaxfockae``: Defaults to None.
         - ``nomega``: Defaults to None.
         - ``precfock``: Defaults to None.
         - ``encutgw``: Defaults to None.
         - ``lrpa``: Defaults to None.
         - ``lpead``: Defaults to None.
         - ``nelm``: Defaults to None. 
         - ``restart``: the return from previous vasp run to use as restart. 

             >> save_this_object = vasp(some parameters) # makes a vasp call.
             >> # make a second vasp using WAVECAR and whatnot from above call
             >> vasp(other parameters, ..., restart = save_this_object) 

             Checks nonscf attribute to correctly perform non-self-consistent
             calculation when required.

         - `relaxation`: sets degrees of freedom to relax. Easier to use
             than isif, nsw, and friends.
         - `set_smearing`: to easily set sigma and ismear.
         - `set_symmetries`: to easily set isym and symprec.

      These parameters can be modified as in C{vasp.ispin = 2} and so forth.
      In the special case that None is given (e.g. C{vasp.ispin = None}), then
      that parameter will not be printed to the INCAR, which means VASP default
      will be used.

      New parameters can be added as follows:

      >>> vasp.add_param = "ispin", 2

      This will print "ISPIN = 2" to the incar. Uppercase conversion is automatic.

      :note: U and NLEP parameters should be given when defining the pseudo-potentials.

      For hackers, how this code works: VASP parameters are set either in
      self.params, in which case VASP key/value pair are the key (in uppercase)
      and value of that dictionary, or in the self.special dictionary. In the
      latter case, the values in the dictionary of classes inheriting from
      `_params.SpecialVaspParam`. They contain (at least) one variable,
      "value", and one method, "incar_string".  The former is the one which is
      modified when an attribute is gotten or set (as in
      C{vasp.some_special_parameter = whatever}). The latter is called, with at
      least vasp as first argument (others are optional, as well as keyword
      args), and should return a valid INCAR string.  
  """ 


  def __init__(self): 
    super(Incar, self).__init__()

    # first, actually sets these two variables by hand, since they are used in __setattr__.
    super(Incar, self).__setattr__("params", {})
    super(Incar, self).__setattr__("special", {})

    self.add_param = "ispin",       1 
    self.add_param = "isif",        0
    self.add_param = "ismear",      None
    self.add_param = "sigma",       None
    self.add_param = "nsw",         None
    self.add_param = "ibrion",      None
    self.add_param = "potim",       None
    self.add_param = "nbands",      None
    self.add_param = "lorbit",      None
    self.add_param = "lsorbit",     None
    self.add_param = "lplane",      None
    self.add_param = "addgrid",     None
    self.add_param = "isym",        None
    self.add_param = "symprec",     None
    self.add_param = "lcorr",       None
    self.add_param = "nupdown",     None
    self.add_param = "lmaxmix",     4
    self.add_param = "weimin",      0
    self.add_param = "lmaxfockae",  None
    self.add_param = "nomega",      None
    self.add_param = "encutgw",     None
    self.add_param = "encutlf",     None
    # objects derived from SpecialVaspParams will be recognized as such and can
    # be added without further fuss.
    self.nelect      = NElect(0)
    self.algo        = Algo("fast")
    self.precision   = Precision("accurate")
    self.ediff       = Ediff(1e-4)
    self.ediffg      = Ediff(None, "ediffg")
    self.encut       = Encut(None)
    self.fftgrid     = FFTGrid(None)
    self.restart     = Restart(None)
    self.U_verbosity = UParams("occupancy")
    self.iniwave     = IniWave(None)
    self.magmom      = Magmom()
    self.npar        = Npar(None)
    self.lwave       = Boolean("lwave", False)
    self.precfock    = PrecFock(None)
    self.lrpa        = Boolean("lrpa", None)
    self.loptics     = Boolean("loptics", None)
    self.lpead       = Boolean("lpead", None)
    self.nelm        = Integer("nelm", None)
    self.nonscf      = NonScf(False)


  def incar_lines(self, *args, **kwargs):
    """ List of incar lines. """

    # prints system name. This is not an option!
    result = []
    if hasattr(self._system, "name"):
      if len(self._system.name) != 0:
        result.append("SYSTEM = \"%s\"\n" % (self._system.name))
    # gathers special parameters.
    # Calls them first in case they change normal key/value pairs.
    specials = []
    for key, value in self.special.items():
      if value.value is None: continue
      line = value.incar_string(self, *args, **kwargs)
      if line is not None: specials.append(line + "\n")
    # prints key/value pairs
    for key, value in self.params.items():
      if value is None: continue
      if isinstance(value, bool):  value = ".TRUE." if value else ".FALSE."
      else: 
        try: value = str(value)
        except ValueError: 
          print "Could not convert vasp parameter %s to string: " % (key), value, "."
          raise
      result.append( "%18s = %s\n" % (key.upper(), value))
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
    elif name in self.params:
      self.params[name] = value
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
  def set_symmetries(self, value):
    """ Type of symmetry used in the calculation.
  
        This sets isym and symprec vasp tags.
        Can be "off" or a float corresponding to the tolerance used to determine
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
  def set_smearing(self, args):
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
      self.ismear, self.sigma = None, None
      return

    if isinstance(args, str): args = args.split()
    first = args[0]
    has_second = len(args) > 1
    second = args[1] if len(args) > 1 else None
    has_third = len(args) > 2
    third = args[2] if len(args) > 2 else None

    if first is None:
      self.ismear = None
      if has_second: self.sigma = None
    elif first == "fermi" or first == "-1":    
      self.ismear = -1
      if has_second: self.sigma = second
    elif first == "gaussian" or first == "0":
      self.ismear = 0
      if has_second: self.sigma = second
    elif first == "metal":
      self.ismear = 1
      if has_second: self.sigma = second
    elif first == "mp" or first == "metal":
      if has_third:
        self.ismear = second
        self.sigma = third
      elif has_second: self.sigma = second
      else: self.ismear = 1
      assert self.ismear >= 1, "Mehtfessel-Paxton order must be at least 1."
    elif first == "bloechl" or first == "-5" or first == "insulator":
      self.ismear = -5
      if has_second: self.sigma = second
    elif first == "tetra" or first == "-4":
      self.ismear = -4
      if has_second: self.sigma = second
    else: 
      try: self._value = int(first)
      except: raise ValueError, "Unknown smearing value %s.\n" % (value)
      assert self._value >= 1, "Unknown smearing value %s.\n" % (value)

  @property
  def relaxation(self):
    """ Sets type of relaxation.
    
        It accepts a tuple, as in:

        >>> vasp.relaxation = "static", 
      
        Some of the parameters (purposefully left out above) are optional:
        
        - first argument can be "static", or a combination of "ionic",
          "cellshape", and "volume".  
        - second (optional) argument is nsw
        - third (optional) argument is ibrion
        - fourth (optional) argument is potim.
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
    return result.lstrip()

  @relaxation.setter
  def relaxation(self, args): 
    """ Returns the kind of relaxation being performed. """
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
      self.params["isif"] = 2
      self.params["ibrion"] = -1
      assert ibrion is None or ibrion == -1, \
             ValueError("Cannot set ibrion to anything but -1 for static calculations.")
      assert nsw  is None or nsw == 0, \
             ValueError("static calculation with nsw > 0 is way too creative.")
      self.params["nsw"] = None
      if potim is not None: self.params["potim"] = potim

    else: # Some kind of relaxations. 
      # ionic calculation.
      if ionic and (not cellshape) and (not volume):   self.params["isif"] = 2
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
        assert ibrion != -1, ValueError("Cannot set ibrion to -1 with strain relaxations.")
        assert ibrion != 0 or self.params["isif"] == 1,\
               ValueError("Cannot set ibrion to 0 with strain relaxations.")
        self.params["ibrion"] = ibrion
      if nsw is not None:
        assert nsw > 0, ValueError("Cannot set nse < 1 and perform strain relaxations.")
        self.params["nsw"] = nsw
      elif self.params["nsw"] is None or self.params["nsw"] == 0: self.params["nsw"] = 50
      if potim is not None: self.params["potim"] = potim
      if self.ediffg is not None:
        if self.ediffg < self.ediff and self.ediffg > 0: self.ediffg = None

  @add_setter
  def set_relaxation(self, value):
    """ Alias to relaxation. """
    from warnings import warn
    warn( DeprecationWarning("set_relaxation is obsolete. Please use relaxation instead."), 
          stacklevel=2 )
    self.relaxation = value

  def __iter__(self):
    """ Iterates over attribute names and values. """
    for key, value in self.params: yield key, value
    for key, value in self.special: yield key, value.value

  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("params")
    special = d.pop("special")
    return d, params, special
  def __setstate__(self, args):
    super(Incar, self).__setattr__("params", args[1])
    super(Incar, self).__setattr__("special", args[2])
    d = self.__dict__.update(args[0])
