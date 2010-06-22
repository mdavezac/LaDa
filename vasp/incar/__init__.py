""" Subpackage defining vasp incar parameters. """
from _params import SpecialVaspParam, NElect, Algo, Precision, Ediff,\
                    Encut, FFTGrid, Restart, UParams
from ...opt.decorators import add_setter

class Incar(object):
  """ Contains vasp Incar parameters. 

      The following assumes you know how to write an INCAR. Although you won't
      need to anymore.  This class separates vasp parameters from methods to
      launch and control vasp. vasp attributes can be listed by iterating over
      this class, or calling iter.

      Parameters(attributes) which are present by default are the following:
         - ispin: Sets number of spins. Must be either 1 or 2. 
         - ismear: Smearing function. Can be set with property L{smearing}. 
         - sigma: Smearing parameter. Can be set with property L{smearing}.
         - isif: Degrees of freedom to relax. Can be set using L{self.relaxation}. 
         - nsw: Number of ionic steps. Can be set using L{self.relaxation}. 
         - ibrion: ionic-relaxation method. Can be set using L{self.relaxation}. 
         - potim: ionic-relaxation step. Can be set using L{self.relaxation}. 
         - L{iniwave}: initial wavefunction to use can be either "random" or "jellium".   
         - nelect: sets number of electrons in calculation above and beyond valence.
             - 0(default) lets VASP compute it from species in the system. 
             - 1 would mean +1 electron
             - -1 would mean +1 hole
             - etc
         - algo: electronic minimization. Can be \"very fast\", \"fast\", or \"normal\" (default). 
         - precision: sets accuracy of calculation. Can be \"accurate\"
           (default), \"low\", \"medium\", \"high\".
         - ediff: sets the tolerance per atom of electronic minimization.
            This tolerance is multiplied by the number of atoms in the system,
            eg consistent from one system to another.
         - encut: factor by which ENMAX of species is multiplied.
         - fftgrid: Either a 3-tuple of integers setting NGX and friend, or
           anything else (other than None). In the latter case, a fake VASP
           calculation is performed to get VASP recommended values.
         - restart: the return from previous vasp run to use as restart. 

             >> save_this_object = vasp(some parameters) # makes a vasp call.
             >> # make a second vasp using WAVECAR and whatnot from above call
             >> vasp(other parameters, ..., restart = save_this_object) 

         - C{set_relaxation}: sets degrees of freedom to relax. Easier to use
             than isif, nsw, and friends.
         - C{set_smearing}: to easily set sigma and ismear.
         - C{set_symmetries}: to easily set isym and symprec.

      These parameters can be modified as in C{vasp.ispin = 2} and so forth.
      In the special case that None is given (e.g. C{vasp.ispin = None}), then
      that parameter will not be printed to the INCAR, which means VASP default
      will be used.

      New parameters can be added as follows:

         >>> vasp.add_param = "ispin", 2

      This will print "ISPIN = 2" to the incar. Uppercase conversion is automatic.

      @note: U and NLEP parameters should be given when defining the pseudo-potentials.

      For hackers, how this code works: VASP parameters are set either in
      self.params, in which case VASP key/value pair are the key (in uppercase)
      and value of that dictionary, or in the self.special dictionary. In the
      latter case, the values in the dictionary of classes inheriting from
      LL{_params.SpecialVaspParam}. They contain (at least) one variable,
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
#    # Then makes sure epydoc reads this right.
#   self.params = {} 
#   """ Key/Value vasp pairs.
#   
#       The vasp key is the uppercase of the dictionary key.
#       In other words C{self.params["ispin"] = 2} will print as "ISPIN = 2".
#   """
#   self.special = {} 
#   """ Special vasp parameters.
#   
#       These parameters need to know about other vasp parameters and/or the
#       system of interest to print INCAR string. 
#   """ 

    self.iniwave = "random"
    self.add_param = "ispin",       1 
    self.add_param = "isif",        0
    self.add_param = "ismear",      None
    self.add_param = "sigma",       None
    self.add_param = "nsw",         None
    self.add_param = "ibrion",      None
    self.add_param = "potim",       None
    self.add_param = "nbands",      None
    self.add_param = "lorbit",      None
    self.add_param = "npar",        None
    self.add_param = "lplane",      None
    self.add_param = "addgrid",     None
    self.add_param = "isym",        None
    self.add_param = "symprec",     None
    self.add_param = "lcorr",       None
    self.add_param = "magmom",      None
    self.add_param = "nelect",      NElect(0)
    self.add_param = "algo",        Algo("normal")
    self.add_param = "precision",   Precision("accurate")
    self.add_param = "ediff",       Ediff(1e-4)
    self.add_param = "encut",       Encut(None)
    self.add_param = "fftgrid",     FFTGrid(None)
    self.add_param = "restart",     Restart(None)
    self.add_param = "U_verbosity", UParams("occupancy")



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
      if value.value == None: continue
      line = value.incar_string(self, *args, **kwargs)
      if line != None: specials.append(line + "\n")
    # prints key/value pairs
    for key, value in self.params.items():
      if value == None: continue
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
    if isinstance(value, SpecialVaspParam): self.special[key] = value
    else: self.params[key] = value

  def __getattr__(self, name): 
    """ Gets a VASP parameter from standard and special dictionaries. """
    if name in self.params: return self.params[name]
    elif name in self.special: return self.special[name].value
    raise AttributeError("Unknown parameter " + name)

  def __setattr__(self, name, value):
    """ Sets a VASP parameter to standard and special dictionaries. """
    if   name in self.params: self.params[name] = value
    elif name in self.special: self.special[name].value = value
    else: super(Incar, self).__setattr__(name, value)

  def __delattr__(self, name): 
    """ Deletes a VASP parameter from standard and special dictionaries. """
    if name in self.params: return self.params.pop(name)
    elif name in self.params: return self.special.pop(name).value
    else: super(Incar, self).__detattr__(name)

  def _get_iniwave(self):
    """ Initializes wave functions with \"random\" or 1(default), \"jellium\" or 2. """ 
    return self.params["iniwave"] 
  def _set_iniwave(self, value):
    value = str(value).split()[0] # removes spaces.
    if value == "1" or value == "random": self.params["iniwave"] = 2
    elif value == "2" or value == "jellium": self.params["iniwave"] = 1
    else: raise ValueError("iniwave cannot be set to " + value + ".")
  iniwave = property(_get_iniwave, _set_iniwave)

  @add_setter
  def set_symmetries(self, value):
    """ Type of symmetry used in the calculation.
  
        This sets isym and symprec vasp tags.
        Can be "off" or a float corresponding to the tolerance used to determine
        symmetry operation. 
    """
    if value == None: self.isym = None
    if str(value).lower() == "off" or str(value) == "0": self.params["isym"] = 0
    elif "isym" in self.params:
      if self.isym == 0: self.isym = None
      self.symprec = value
    else: self.symprec = value

  @add_setter
  def set_smearing(self, args):
    """ Value of the smearing used in the calculation. 
  
        It can be specified as a string L{vasp.smearing = "type", x}, where type
        is any of fermi, gaussian, mp, tetra, metal or insulator, and x is the
        energy scale in eV.
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
    if args == None: 
      self.ismear, self.sigma = None, None
      return

    first = args[0]
    has_second = len(args) > 1
    second = args[1] if len(args) > 1 else None
    has_third = len(args) > 2
    third = args[2] if len(args) > 2 else None

    if first == None:
      self.ismear = None
      if has_second: self.sigma = None
    elif first == "fermi" or first == "-1":    
      self.ismear == -1
      if has_second: self.sigma = second
    elif first == "gaussian" or first == "0":
      self.ismear == 0
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

  @add_setter
  def set_relaxation(self, *args): 
    """ Sets type of relaxation.
    
          - first argument can be "static", or a combination of "ion(ic|s)",
              "cell(\s+|-|_?(?:shape)?", and "volume". 
            Can also be set using an integer between 0 and 7. See VASP manual. 
          - second (optional) argument is nsw
          - third (optional) argument is ibrion
          - fourth (optional) argument is potim.
    """
    import re
    if args[0] == None: isif = None
    else:
      try: isif = int(args[0]) 
      except ValueError: 
        isif = str(args[0]).lower()
        ionic = re.search( "ion(ic|s)?", isif.lower() ) != None
        cellshape = re.search( "cell(\s+|-|_)?(?:shape)?", isif.lower() ) != None
        volume = re.search( "volume", isif.lower() ) != None
        if (not ionic) and (not cellshape) and (not volume): isif = 0
        elif ionic and (not cellshape) and (not volume):     isif = 2
        elif ionic and cellshape and volume:                 isif = 3
        elif ionic and cellshape and (not volume):           isif = 4
        elif(not ionic) and  cellshape and (not volume):     isif = 5
        elif(not ionic) and  cellshape and volume:           isif = 6
        elif(not ionic) and (not cellshape) and volume:      isif = 7
        elif ionic and (not cellshape) and volume:
          raise RuntimeError, "VASP does not allow relaxation of atomic position"\
                              "and volume at constant cell-shape.\n"
    self.params["isif"] = isif
    if len(args) > 1: self.params["nsw"] = args[1]
    if len(args) > 2: self.params["ibrion"] = args[2]
    if len(args) > 3: self.params["potim"] = args[3]

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
