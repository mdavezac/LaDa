""" Standard parameter types for use as attributes in Incar """
__docformat__ = "restructuredtext en"
__all__ = [ "SpecialVaspParam", "NElect", "Algo", "Precision", "Ediff",\
            "Encut", "FFTGrid", "Restart", "UParams", "IniWave",\
            "Magmom", 'Npar', 'Boolean', 'Integer', 'Choices', 'PrecFock', 'NonScf']
class SpecialVaspParam(object): 
  """ Type checking class. """
  def __init__(self, value): 
    super(SpecialVaspParam, self).__init__()
    self.value = value
  def __repr__(self): return "{0.__class__.__name__}({1})".format(self, repr(self.value))

class Magmom(SpecialVaspParam):
  """ Prints magmom to INCAR. 

      There are three types of usage: 
        - do nothing if the instance's value is None or False.
        - print a string preceded by "MAGMOM = " if the instance's value is a string. 
        - print the actual MAGMOM string from the magnetic moments attributes
          ``magmom`` in the structure's atoms.
  """
  def __init__(self, value = "attribute: magmom"):
    """ Initializes magmom. """
    super(Magmom, self).__init__(value)
    
  def incar_string(self, *args, **kwargs):
    """ Prints the magmom string if requested. """
    from ...crystal import specieset
    if self.value is None or self.value == False: return None
    if kwargs["vasp"].ispin == 1: return None
    if isinstance(self.value, str): return "MAGMOM = {0}".format(self.value)
    
    structure = kwargs['structure']
    result = ""
    for specie in specieset(structure):
      moments = [getattr(u, 'magmom', 0e0) for u in structure if u.type == specie]
      tupled = [[1, moments[0]]]
      for m in moments[1:]: 
        if abs(m - tupled[-1][1]) < 1e-12: tupled[-1][0] += 1
        else: tupled.append([1, m])
      for i, m in tupled:
        if i == 1: result += "{0:.2f} ".format(m)
        else:      result += "{0}*{1:.2f} ".format(i, m)
    return 'MAGMOM = {0}'.format(result)
  
class Npar(SpecialVaspParam):
  """ Parallelization over bands. 

      This parameter is described `here`__.
      It can be set to a particular number:
  
      >>> vasp.npar = 2

      Or it can be deduced automatically. In the latter case, npar is set to the
      largest power of 2 which divides the number of processors:
 
      >>> vasp.npar = "power of two"
      
      If the number of processors is not a power of two, prints nothing.

      .. __: http://cms.mpi.univie.ac.at/vasp/guide/node138.html
  """

  def __init__(self, value): super(Npar, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    from re import search
    from math import log, sqrt
    if self.value is None: return None
    if not isinstance(self.value, str): 
      if self.value < 1: return None
      return "NPAR = {0}".format(self.value)

    comm = kwargs.get('comm', None)
    n = getattr(comm, 'n', getattr(comm, 'size', -1))
    if n == 1: return None
    if self.value == "power of two":
      m = int(log(n)/log(2))
      for i in range(m, 0, -1):
        if n % 2**i == 0: return "NPAR = {0}".format(i)
      return None
    if self.value == "sqrt": return "NPAR = {0}".format(int(sqrt(n)+0.001))
    raise ValueError("Unknown request npar = {0}".format(self.value))
    
class NElect(SpecialVaspParam):
  """ Sets number of electrons relative to neutral system.
      
      Gets the number of electrons in the (neutral) system. Then adds value to
      it and computes with the resulting number of electrons.
      >>> nelect = NElect(0) # charge neutral system
      >>> nelect.value = 1   # charge -1 (1 extra electron)
      >>> nelect.value = -1  # charge +1 (1 extra hole)

      :Param value: (default:0) number of electrons to add to charge neutral
                    system.
  """

  def __init__(self, value): super(NElect, self).__init__(value)

  def nelectrons(self, vasp, structure):
    """ Total number of electrons in the system """
    from math import fsum
    # constructs dictionnary of valence charge
    valence = {}
    for key, value in vasp.species.items():
      valence[key] = value.valence
    # sums up charge.
    return fsum( valence[atom.type] for atom in structure )
    
  def incar_string(self, *args, **kwargs):
    # gets number of electrons.
    charge_neutral = self.nelectrons(kwargs['vasp'], kwargs['structure'])
    # then prints incar string.
    if self.value == 0:
      return "# NELECT = {0} Charge neutral system".format(charge_neutral)
    elif self.value > 0:
      return "NELECT = {0}  # negatively charged system ({1})"\
             .format(charge_neutral + self.value, -self.value)
    else: 
      return "NELECT = {0}  # positively charged system (+{1})"\
             .format (charge_neutral + self.value, -self.value)
          
      
class Algo(SpecialVaspParam): 
  """ Electronic minimization. 
  
      Defines the `ALGO`__ tag.
      Takes one of the following values:
        - very fast
        - fast, f (default)
        - normal, n
        - all, a
        - damped, d 
        - Diag 
        - conjugate, c (vasp 5)
        - subrot (vasp 5)
        - eigenval (vasp 5)
        - Nothing (vasp 5)
        - Exact  (vasp 5)
        - chi
        - gw
        - gw0
        - scgw
        - scgw0

      .. warning:: The string None is not  allowed, as it would lead to
         confusion with the python object None. Please use "Nothing" instead.
         The python object None will simply not print the ALGO keyword to the
         INCAR file.

      .. note:: By special request, "fast" is the default algorithm.

      .. __: http://cms.mpi.univie.ac.at/vasp/vasp/ALGO_tag.html
  """ 
  def __init__(self, value="fast"): super(Algo, self).__init__(value)
  @property
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return None
    try: from lada import is_vasp_4
    except: is_vasp_4 == False
    if not hasattr(value, 'lower'):
      raise TypeError("ALGO cannot be set with {0}.".format(value))
    lower = value.lower().rstrip().lstrip()
    lower = lower.replace('_', '')
    lower = lower.replace('-', '')
    if is_vasp_4                                         \
       and ( lower[0] in ['c', 's', 'e']                 \
             or lower in [ "nothing", "subrot", "exact", \
                           "gw", "gw0", "chi", "scgw",   \
                           "scgw0"] ): 
      raise ValueError, "algo value ({0}) is not valid with VASP 4.\n".format(value)

    if lower == "diag": value = "Diag"
    elif lower == "nothing": value = "Nothing"
    elif lower == "chi":  value = "chi"
    elif lower == "gw":   value = "GW"
    elif lower == "gw0":  value = "GW0"
    elif lower == "scgw": value = "scGW"
    elif lower == "scgw0": value = "scGW0"
    elif lower[0] == 'v': value = "Very_Fast" if is_vasp_4 else 'VeryFast'
    elif lower[0] == 'f': value = "Fast"
    elif lower[0] == 'n': value = "Normal"
    elif lower[0] == 'd': value = "Damped"
    elif lower[0] == 'a': value = "All"
    elif lower[0] == 'c': value = "Conjugate"
    elif lower[0] == 's': value = "Subrot"
    elif lower[0] == 'e': value = "Eigenval"
    else:
      self._value = None
      raise ValueError("algo value ({0!r}) is invalid.\n".format(value))
    self._value = value
    
  def incar_string(self, *args, **kwargs):
    if self.value is None: return None
    return "ALGO = {0}".format(self.value)

class Precision(SpecialVaspParam):
  """ Sets accuracy of calculation. 

      Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
  """
  def __init__(self, value): super(Precision, self).__init__(value)
  def incar_string(self, *args, **kwargs):
    value = self.value.lower()
    if value not in ["accurate", "lower", "medium", "high"]:
      raise ValueError("PRECISION value ({0}) is not allowed.".format(self.key))
    return "PREC = {0}".format(self.value)

class Ediff(SpecialVaspParam):
  """ Sets the convergence criteria (per atom) for electronic minimization.

      This tolerance is multiplied by the number of atoms in the system. This
      makes tolerance consistent from one system to the next.
  """
  def __init__(self, value, name="ediff"):
    """ Creates *per atom* tolerance. """
    super(Ediff, self).__init__(value)
    self.name = name
  def incar_string(self, *args, **kwargs):
    if self.value is None: return 
    if self.value < 0: 
      return "{0} = {1} ".format(getattr(self, "name", "ediff").upper(), self.value)
    return "{0} = {1} ".format( getattr(self, "name", "ediff").upper(),
                                self.value * float(len(kwargs["structure"])) )
  def __repr__(self):
    """ Representation of Ediff. """
    name = getattr(self, "name", "ediff")
    if name == "ediff": return "{0.__class__.__name__}({1})".format(self, repr(self.value))
    return "{0.__class__.__name__}({1}, {2})".format(self, repr(self.value), repr(name))

class Encut(SpecialVaspParam):
  """ Defines cutoff factor for calculation. 

      There are three ways to set this parameter:

      - if 0 < value <= 3, then the cutoff is value * ENMAX, where ENMAX is
        the maximum recommended cutoff for the species in the system.
      - if value > 3, then prints encut is exactly value.
      - if 0 or None, does not print anything to INCAR

      If 0 or None, uses VASP default.
  """
  def __init__(self, value): super(Encut, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    """ Prints ENCUT parameter. """
    from math import fabs
    from ...crystal import specieset
    assert self.value > -1e-12, ValueError("Wrong value for cutoff.")
    if fabs(self.value) < 1e-12: return "# ENCUT = VASP default"
    if self.value > 3e0 + 1e-12:
      return "ENCUT = {0}".format(self.value.rescale("eV") if hasattr(self.value, "rescale") else self.value)
    types = specieset(kwargs["structure"])
    encut = max(kwargs["vasp"].species[type].enmax for type in types)
    return "ENCUT = {0} ".format(float(encut) * self.value)

class FFTGrid(SpecialVaspParam):
  """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """
  def __init__(self, value): super(FFTGrid, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    if self.value is None: return None
    return "NGX = {0[0]}\nNGY = {0[1]}\nNGZ = {0[2]}".format(self.value)

class Restart(SpecialVaspParam):
  """ Return from previous run from which to restart.
      
      If None, then starts from scratch.
  """
  def __init__(self, value): super(Restart, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    from os.path import join, exists, getsize
    from shutil import copy
    from ...misc import copyfile
    from .. import files
    nonscf = getattr(kwargs["vasp"], 'nonscf', False)
    istart = "0   # start from scratch"
    icharg = "{0}   # superpositions of atomic densities".format(12 if nonscf else 2)
    if self.value is None: istart = "0   # start from scratch"
    elif not self.value.success:
      istart = "0   # start from scratch"
    else:
      ewave = exists( join(self.value.directory, files.WAVECAR) )
      if ewave: ewave = getsize(join(self.value.directory, files.WAVECAR)) > 0
      echarge = exists( join(self.value.directory, files.CHGCAR) )
      if echarge: echarge = getsize(join(self.value.directory, files.CHGCAR)) > 0
      if ewave:
        path = join(self.value.directory, files.WAVECAR)
        istart = "1   # restart"
        icharg = "{0}   # from wavefunctions {1}".format(10 if nonscf else 0, path)
        copy(path, ".")
        if echarge: copy(join(self.value.directory, files.CHGCAR), '.')
      elif echarge:
        path = join(self.value.directory, files.CHGCAR)
        istart = "1   # restart"
        icharg = "{0}   # from charge {1}".format(11 if nonscf else 1, path)
        copy(path, ".")
      else: 
        istart = "0   # start from scratch"
        icharg = "{0}   # superpositions of atomic densities".format(12 if nonscf else 2)
      copyfile(join(self.value.directory, files.EIGENVALUES), nothrow='same exists',
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.CONTCAR), files.POSCAR,\
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),\
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.WAVEDER), files.WAVEDER,
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),
               nocopyempty=True) 
    return  "ISTART = {0}\nICHARG = {1}".format(istart, icharg)

class NonScf(SpecialVaspParam):
  """ Return from previous run from which to restart.
      
      If None, then starts from scratch.
  """
  def __init__(self, value):  super(NonScf, self).__init__(value)
  @property
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if isinstance(value, str):
      if len(value) == 0: value = False
      elif   value.lower() == "true"[:min(len(value), len("true"))]: value = True
      elif value.lower() == "false"[:min(len(value), len("false"))]: value = False
      else: raise RuntimeError("Uknown value for nonscf: {0}").format(value)
    self._value = value == True

  def incar_string(self, *args, **kwargs): return None
  def __repr__(self): return "{0.__class__.__name__}({0.value!r})".format(self)

class UParams(SpecialVaspParam): 
  """ Prints U parameters if any found in species settings """
  def __init__(self, value):
    import re
    
    if value is None: value = 0
    elif hasattr(value, "lower"): 
      value = value.lower() 
      if value == "off": value = 0
      elif value == "on": value = 1
      elif None != re.match(r"\s*occ(upancy)?\s*", value): value = 1
      elif None != re.match(r"\s*(all|pot(ential)?)\s*", value): value = 2

    super(UParams, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    """ Prints LDA+U INCAR parameters. """
    from ...crystal import specieset
    types = specieset(kwargs['structure'])
    species = kwargs['vasp'].species
    # existence and sanity check
    has_U, which_type = False, None 
    for type in types:
      specie = species[type]
      if len(specie.U) == 0: continue
      if len(specie.U) > 4: 
        raise AssertionError, "More than 4 channels for U/NLEP parameters"
      has_U = True
      # checks consistency.
      which_type = specie.U[0]["type"]
      for l in specie.U[1:]: 
        assert which_type == l["type"], \
               AssertionError("LDA+U/NLEP types are not consistent across species.")
    if not has_U: return "# no LDA+U/NLEP parameters";

    # Prints LDA + U parameters
    result = "LDAU = .TRUE.\nLDAUPRINT = {0}\nLDAUTYPE = {1}\n".format(self.value, which_type)

    for i in range( max(len(species[type].U) for type in types) ):
      line = "LDUL{0}=".format(i+1), "LDUU{0}=".format(i+1), "LDUJ{0}=".format(i+1), "LDUO{0}=".format(i+1)
      for type in types:
        specie = species[type]
        a = -1, 0e0, 0e0, 1
        if len(specie.U) <= i: pass
        elif specie.U[i]["func"] == "U":    
          a = [specie.U[i]["l"], specie.U[i]["U"], specie.U[i]["J"], 1]
        elif specie.U[i]["func"] == "nlep": 
          a = [specie.U[i]["l"], specie.U[i]["U0"], 0e0, 2]
        elif specie.U[i]["func"] == "enlep":
          a = [specie.U[i]["l"], specie.U[i]["U0"], specie.U[i]["U1"], 3]
        else: raise RuntimeError, "Debug Error."
        if hasattr(a[1], "rescale"): a[1] = a[1].rescale("eV")
        if hasattr(a[2], "rescale"): a[2] = a[2].rescale("eV")
        line = "{0[0]} {1[0]}".        format(line, a),\
               "{0[1]} {1[1]:18.10e}". format(line, a),\
               "{0[2]} {1[2]:18.10e}".format(line, a),\
               "{0[3]} {1[3]}".        format(line, a)
      result += "\n{0}\n{1}\n{2}\n{3}\n".format(*line)
    return result
  def __repr__(self):
    """ Representation of UParams """
    return "{0.__class__.__name__}({1!r})".format(self, ["off", "on", "all"][self.value])

class IniWave(SpecialVaspParam):
  def __init__(self, value): super(IniWave, self).__init__(value)
  def incar_string(self, *args, **kwargs):
    """ Returns VASP incar string. """
    if self.value == "1" or self.value == "random": result = 1
    elif self.value == "0" or self.value == "jellium": result = 0
    else: raise ValueError("iniwave cannot be set to " + self.value + ".")
    return "INIWAVE = {0}\n".format(result)


class Boolean(SpecialVaspParam):
  """ Any boolean vasp parameters. 
  
      Python is very liberal in how it converts any object to a boolean, eg an
      empty dictionary is false while non-empty dictionary is true.
      In order to keep this behavior, the value given to this parameter is kept
      as is as long as possible, and converted only when writing the incar. The
      only difference with the python behavior is that if using strings (which generally
      evaluate to true or depending whether or not they are empty), these must
      be "True" or "False", or variations thereoff. The empty string will
      evaluate to the VASP default (eg equivalent to using None).
  """
  def __init__(self, key, value):
    super(Boolean, self).__init__(value)
    self.key = key
    """ VASP key corresponding to this input. """
  @property
  def value(self):
    return self._value
  @value.setter
  def value(self, value):
    if isinstance(value, str):
      if len(value) == 0: value = False
      elif   value.lower() == "true"[:min(len(value), len("true"))]: value = True
      elif value.lower() == "false"[:min(len(value), len("false"))]: value = False
      else: raise TypeError("Cannot interpret string {0} as a boolean.".format(value))
    self._value = value == True
  def incar_string(self, *args, **kwargs):
    value = self._value
    if isinstance(value, str):
      if len(value) == 0: value is None 
      elif value.lower() == "true"[:len(value)]: value = True
      else: value = False
    if self.value is None: return None
    return "{0} = {1}".format(self.key.upper(), ".TRUE." if bool(self.value) else ".FALSE.")
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1!r}, {2!r})".format(self, self.key, self.value)

class Integer(SpecialVaspParam):
  """ Any integer vasp parameters. 
  
      The value is always of type integer. Other types are converted to an
      integer where possible, and will throw TypeError otherwise.
  """
  def __init__(self, key, value):
    super(Integer, self).__init__(value)
    self.key = key
    """ VASP key corresponding to this input. """
  @property 
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return
    try: self._value = int(value)
    except: raise TypeError("Could not evaluate {0} as an integer.".format(value))
  def incar_string(self, *args, **kwargs):
    if self.value is None: return None
    return "{0} = {1}".format(self.key.upper(), self.value)
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1}, {2})".format(self, repr(self.key), repr(self.value))

class Choices(SpecialVaspParam):
  """ Vasp parameters with a limited set of choices. """
  def __init__(self, key, choices, default=None):
    """ Initializes the Choices-type vasp parameters.

        :param key:
	    Name of the VASP parameter, e.g. "precfock". It needs not be in
            uppercase. In fact, lower case is preferred for being more pythonic.
        :param choices:
            Dictionary where key is an allowed VASP input for this parameter.
            To each key is associated a list (or set), with allowable forms
            which will translate to the key in the incar. A modified copy of
            this dictionary is owned by the instance being initialized. All
            keys and items should be meaningfully convertible to strings.
        :param default :
            Option from ``choices`` to use as default.

      .. note:: The keys are case-sensitive. The values are not.
    """
    self.key = key
    """ VASP key corresponding to this input. """
    self.choices = {}
    """ Allowable set of choices. """
    for key, items in choices.iteritems():
      self.choices[key] = set( [u.lower() if hasattr(u, 'lower') else u for u in items]\
                               + [key, str(key).lower()])
    super(Choices, self).__init__(default)

  @property
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return
    if hasattr(value, 'lower'): value == value.lower()
    for key, items in self.choices.iteritems():
      if value in items: self._value = key; return
    raise ValueError("{0} is not an acceptable choice for {1.key}: {1.choices}.".format(value, self))
  def incar_string(self, *args, **kwargs):
    if self.value is None: return None
    return "{0} = {1}".format(self.key.upper(), self.value)
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1}, {2}, {3})"\
           .format(self, repr(self.key), repr(self.choices), repr(self.value))

class PrecFock(Choices):
  """ Sets up PRECFOCK parameter. 
      
      Allowable options are:

      - L or low:      coarse grid for HF, normal augmentation charge.
      - M or medium:   normal grid for HF, normal augmentation charge.
      - F or fast:     coarse grid for HF, soft augmentation charge. 
      - N or normal:   PREC=N grid for HF, soft augmentation charge. 
      - A or accurate: PREC=A grid for HF, soft augmentation charge.

      .. note:: The values are not case-sensitive. 
  """
  def __init__(self, value=None):
    """ Initializes PRECFOCK parameter. """
    choices = { 'L': ['low'], 'M': ['medium'], 'F': ['fast'],
                'N': ['normal'], 'A': ['accurate'] }
    super(PrecFock, self).__init__("precfock", choices, value)
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1})".format(self, repr(self.value))
