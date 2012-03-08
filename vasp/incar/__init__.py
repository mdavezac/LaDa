""" Subpackage defining vasp incar parameters. """
__docformat__ = "restructuredtext en"
__all__ = [ "SpecialVaspParam", "NElect", "Algo", "Precision", "Ediff",\
            "Encut", "FFTGrid", "Restart", "UParams", "IniWave", 'Ediffg', "EncutGW", \
            "Incar", "Magmom", 'Npar', 'Boolean', 'Integer', 'Choices', 'PrecFock',
            "System", 'PartialRestart', 'Relaxation' ]
from _params import SpecialVaspParam, NElect, Algo, Precision, Ediff,\
                    Encut, FFTGrid, PartialRestart, Restart, UParams, IniWave, Magmom,\
                    Npar, Boolean, Integer, PrecFock, NonScf, Ediffg, Choices, \
                    EncutGW, System, Relaxation
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
    self.add_param = "ismear",      None
    self.add_param = "isigma",      None
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
    self.system      = System(True)
    self.Relaxation  = Relaxation(None)

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

    # gathers special parameters.
    # Calls them first in case they change normal key/value pairs.
    result, specials, comments = [], [], []
    for key, value in self.special.items():
      if value.value is None: continue
      line = value.incar_string(**kwargs)
      if line is None: continue
      line = line.rstrip().lstrip()
      if line[-1] != '\n': line += '\n'
      if line[0] == '#': comments.append(line); continue
      if '=' in line and line.find('=') < 18:
        line = "{0: <{1}}".format(' ', 19 - line.find('=')) + line
      specials.append(line)
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
    result = sorted(result, key=lambda a: a.lstrip()[0])
    result.extend(comments)
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

  @property
  def symmetries(self):
    """ Type of symmetry used in the calculation.
  
        This sets :py:attr:`isym` and :py:attr:`symprec` vasp tags. Can be
        "off" or a float corresponding to the tolerance used to determine
        symmetry operation. 
    """
    if self.isym is None and self.symprec is None: return True
    if self.isym is None: return self.symprec
    if self.isym == 0: return False

  @symmetries.setter
  def symmetries(self, value):
    if value is None: self.isym = None
    elif str(value).lower() == "off" or value is "0" or value is False: self.params["isym"] = 0
    elif str(value).lower() == "on" or value is True or value is True:
       self.symprec = None
       self.isym = None
    elif isinstance(value, float): 
       self.symprec = value
       self.isym = None
    else: raise ValueError("Uknown value when setting symmetries ({0}).".format(value))

  @property 
  def smearing(self):
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
    ismear = { -1: 'fermi', 0: 'gaussian', 1: 'metal', -5:'bloechl',
               -4: 'tetra', 2: 'mp 2', 3: 'mp 3'}[self.ismear]
    return ismear, self.isigma

  @smearing.setter
  def smearing(self, args):
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
      try: self.ismear = int(first)
      except: raise ValueError("Unknown smearing value {0}.\n".format(first))
      if self.ismear < 1:
        raise ValueError("Unknown smearing value {0}.\n".format(first))


  def __getstate__(self):
    d = self.__dict__.copy()
    params = d.pop("params")
    special = d.pop("special")
    return d, params, special
  def __setstate__(self, args):
    super(Incar, self).__setattr__("params", args[1])
    super(Incar, self).__setattr__("special", args[2])
    d = self.__dict__.update(args[0])
