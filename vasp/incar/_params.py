""" Standard parameter types for use as attributes in Incar """
__docformat__ = "restructuredtext en"
__all__ = [ "SpecialVaspParam", "NElect", "Algo", "Precision", "Ediff",\
            "Ediffg", "Encut", "EncutGW", "EncutLF", "FFTGrid", "Restart", "UParams", "IniWave",\
            "Magmom", 'Npar', 'Boolean', 'Integer', 'Choices', 'PrecFock', 'NonScf']
class SpecialVaspParam(object): 
  """ Base type for special vasp parameters. 
  
      Special vasp parameters do something more than just print to the incar.
      What *more* means depends upon the parameter.

  """
  def __init__(self, value): 
    super(SpecialVaspParam, self).__init__()
    self.value = value
    """ Value derived classes will do something with. """
  def __repr__(self): return "{0.__class__.__name__}({1})".format(self, repr(self.value))

class Magmom(SpecialVaspParam):
  """ Sets the initial magnetic moments on each atom.

      There are three types of usage: 
        - do nothing if the instance's value is None or False.
        - print a string preceded by "MAGMOM = " if the instance's value is a string. 
        - print the actual MAGMOM string from the magnetic moments attributes
          ``magmom`` in the structure's atoms if anything but a string, None,
          or False.

      If the calculation is **not** spin-polarized, then the magnetic moment
      tag is not set.

      .. seealso:: `MAGMOM <http://cms.mpi.univie.ac.at/vasp/guide/node100.html>`_
  """
  def __init__(self, value = "attribute: magmom"):
    super(Magmom, self).__init__(value)
    
  def incar_string(self, **kwargs):
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

      Npar defines how many nodes work on one band.
      It can be set to a particular number:
  
      >>> vasp.npar = 2

      Or it can be deduced automatically. Different schemes are available:
      
        - power of two: npar is set to the largest power of 2 which divides the
          number of processors.
 
          >>> vasp.npar = "power of two"

          If the number of processors is not a power of two, prints nothing.

        - square root: npar is set to the square root of the number of processors.

          >>> vasp.npar = "sqrt"
      

      .. seealso: `NPAR <http://cms.mpi.univie.ac.at/vasp/guide/node138.html>`_
  """

  def __init__(self, value): super(Npar, self).__init__(value)

  def incar_string(self, **kwargs):
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

      :param integer value:
        Number of electrons to add to charge neutral system. Defaults to 0.

      .. seealso:: `NELECT <http://cms.mpi.univie.ac.at/vasp/vasp/NELECT.html>`_
  """

  def __init__(self, value=0): super(NElect, self).__init__(value)

  def nelectrons(self, vasp, structure):
    """ Total number of electrons in the system """
    from math import fsum
    # constructs dictionnary of valence charge
    valence = {}
    for key, value in vasp.species.items():
      valence[key] = value.valence
    # sums up charge.
    return fsum( valence[atom.type] for atom in structure )
    
  def incar_string(self, **kwargs):
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
  
      Defines the kind of algorithm vasp will run.
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

      If :py:data:`is_vasp_4 <lada.is_vasp_4>` is an existing configuration
      variable of :py:mod:`lada` the parameters marked as vasp 5 will fail.

      .. warning:: The string None is not  allowed, as it would lead to
         confusion with the python object None. Please use "Nothing" instead.
         The python object None will simply not print the ALGO keyword to the
         INCAR file.

      .. note:: By special request, "fast" is the default algorithm.

      .. seealso:: `ALGO <http://cms.mpi.univie.ac.at/vasp/vasp/ALGO_tag.html>`_
  """ 
  def __init__(self, value="fast"): super(Algo, self).__init__(value)
  @property
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return None
    try: from lada import is_vasp_4
    except: is_vasp_4 = False
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
    
  def incar_string(self, **kwargs):
    if self.value is None: return None
    return "ALGO = {0}".format(self.value)

class Ediff(SpecialVaspParam):
  """ Sets the convergence criteria (per atom) for electronic minimization.

      - value > 0e0: the tolerance is multiplied by the number of atoms in the
        system. This makes tolerance consistent from one system to the next.
      - value < 0e0: tolerance is given as absolute value, without multiplying
        by size of system.

      .. seealso:: `EDIFF <http://cms.mpi.univie.ac.at/vasp/guide/node105.html>`_
  """
  def __init__(self, value, name="ediff"):
    """ Creates *per atom* tolerance. """
    super(Ediff, self).__init__(value)
    self.name = name
  def incar_string(self, **kwargs):
    if self.value is None: return 
    if self.value < 0: 
      return "{0} = {1} ".format(getattr(self, "name", "ediff").upper(), -self.value)
    return "{0} = {1} ".format( getattr(self, "name", "ediff").upper(),
                                self.value * float(len(kwargs["structure"])) )
  def __repr__(self):
    return "{0.__class__.__name__}({0.value!r})".format(self)

class Ediffg(SpecialVaspParam):
  """ Sets the convergence criteria (per atom) for ionic minimization.

      - value > 0e0: the tolerance is multiplied by the number of atoms in the
        system. This makes tolerance consistent from one system to the next.
      - value < 0e0: tolerance is given as is (negative), and applies to forces.

      .. seealso:: `EDIFFG <http://cms.mpi.univie.ac.at/vasp/guide/node107.html>`_
  """
  def __init__(self, value):
    """ Creates *per atom* tolerance. """
    super(Ediffg, self).__init__(value)
  def incar_string(self, **kwargs):
    if self.value is None: return 
    if self.value < 0: return "EDIFFG = {0} ".format(self.value)
    return "EDIFFG = {0} ".format(self.value * float(len(kwargs["structure"])))
  def __repr__(self):
    return "{0.__class__.__name__}({0.value!r})".format(self)

class Encut(SpecialVaspParam):
  """ Defines cutoff factor for calculation. 

      There are three ways to set this parameter:

      - if value is floating point and 0 < value <= 3: then the cutoff is
        ``value * ENMAX``, where ENMAX is the maximum recommended cutoff for
        the species in the system.
      - if value > 3 eV, then prints encut is exactly value (in eV). Any energy
        unit is acceptable.
      - if value < 0 eV or None, does not print anything to INCAR. 
      
      .. seealso:: `ENCUT <http://cms.mpi.univie.ac.at/vasp/vasp/ENCUT_tag.html>`_
  """
  KEY = "ENCUT"
  def __init__(self, value): super(Encut, self).__init__(value)

  def incar_string(self, **kwargs):
    from math import fabs
    from ...crystal import specieset
    from quantities import eV
    value = self.value
    if value is None: return None
    elif hasattr(self.value, 'rescale'): value = float(value.rescale(eV))
    elif value >= 1e-12 and value <= 3.0:
      types = specieset(kwargs["structure"])
      encut = max(kwargs["vasp"].species[type].enmax for type in types)
      if hasattr(encut, 'rescale'): encut = float(encut.rescale(eV))
      return "{0} = {1} ".format(self.KEY, encut * value)
    if value < 1e-12: return None
    return "{0} = {1}".format(self.KEY, value)

class EncutGW(Encut):
  """ Defines cutoff factor for GW calculation. 

      There are three ways to set this parameter:

      - if value is floating point and 0 < value <= 3: then the cutoff is
        ``value * ENMAX``, where ENMAX is the maximum recommended cutoff for
        the species in the system.
      - if value > 3 eV, then prints encut is exactly value (in eV). Any energy
        unit is acceptable.
      - if value < 0 eV or None, does not print anything to INCAR. 
      
      .. seealso:: `ENCUTGW
        <http://cms.mpi.univie.ac.at/vasp/vasp/ENCUTGW_energy_cutoff_response_function.html>`_
  """
  KEY = "ENCUTGW"
  def __init__(self, value): super(EncutGW, self).__init__(value)

class FFTGrid(SpecialVaspParam):
  """ FFT mesh of the wavefunctions.
  
      This must a sequence of three integers.

      .. seealso:: `NGX, NGY, NGZ
        <http://cms.mpi.univie.ac.at/vasp/guide/node93.html>`_
  """
  def __init__(self, value): super(FFTGrid, self).__init__(value)
  @property 
  def value(self): return self._value
  @value.setter
  def value(self, value): 
    from numpy import array
    if value is None:
      self._value = None
      return
    if len(list(value)) != 3: raise TypeError("FFTGrid expects three numbers.")
    self._value = array(value)
  def incar_string(self, **kwargs):
    if self.value is None: return None
    return "NGX = {0[0]}\nNGY = {0[1]}\nNGZ = {0[2]}".format(self.value)

class Restart(SpecialVaspParam):
  """ Return from previous run from which to restart.
      
      It is either an vasp extraction object of some kind, or None.
      In the latter case, the calculation starts from scratch. 
      However, if an extraction object exists *and* the calculation it refers
      to was successfull, then it will check whether WAVECAR and CHGCAR exist
      and set :py:attr:'istart <incar.Incar.istart>` and :py:attr:`icharg
      <incar.Incar.icharg>` accordingly. It also checks whether
      :py:attr:`nonscf <incar.Incar.nonscf>` is True or False, and sets
      :py:attr:`icharg <incar.Incar.icharg>` accordingly. 

      .. seealso:: `ICHARG
        <http://cms.mpi.univie.ac.at/vasp/guide/node102.html>`_, `ISTART
        <http://cms.mpi.univie.ac.at/vasp/guide/node101.html>`_
  """
  def __init__(self, value): super(Restart, self).__init__(value)

  def incar_string(self, **kwargs):
    from os.path import join, exists, getsize
    from shutil import copy
    from ...misc import copyfile
    from .. import files

    if self.value is None or self.value.success == False:
      if kwargs['vasp'].nonscf: kwargs['vasp'].icharg = 12
      return None
    else:
      ewave = exists( join(self.value.directory, files.WAVECAR) )
      if ewave: ewave = getsize(join(self.value.directory, files.WAVECAR)) > 0
      if ewave:
        copy(join(self.value.directory, files.WAVECAR), ".")
        kwargs['vasp'].istart = 1
      else: kwargs['vasp'].istart = 0
      echarg = exists( join(self.value.directory, files.CHGCAR) )
      if echarg: echarg = getsize(join(self.value.directory, files.CHGCAR)) > 0
      if echarg:
        copy(join(self.value.directory, files.CHGCAR), ".")
        kwargs['vasp'].icharg = 1
      else: kwargs['vasp'].icharg = 0 if kwargs['vasp'].istart == 1 else 2
      if getattr(kwargs["vasp"], 'nonscf', False): kwargs['vasp'].icharg += 10

      copyfile(join(self.value.directory, files.EIGENVALUES), nothrow='same exists',
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.CONTCAR), files.POSCAR,\
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),\
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.WAVEDER), files.WAVEDER,
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.TMPCAR), files.TMPCAR,
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),
               nocopyempty=True) 
    return None

class NonScf(SpecialVaspParam):
  """ Whether to perform a self-consistent or non-self-consistent run. 
  
      Accepts only True or False(default). This parameter works with
      :py:class:`Restart` to determine the value to give :py:attr:`icharg
      <lada.vasp.incar.Incar.icharg>`
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

  def incar_string(self, **kwargs): return None
  def __repr__(self): return "{0.__class__.__name__}({0.value!r})".format(self)

class UParams(SpecialVaspParam): 
  """ Sets U, nlep, and enlep parameters. 
 
      The U, nlep, and enlep parameters of the atomic species are set at the
      same time as the pseudo-potentials. This object merely sets up the incar
      with right input.

      However, it does accept one parameter, which can be "off", "on", "occ" or
      "all" wich defines the level of verbosity of VASP (with respect to the
      parameters).


      .. seealso:: `LDAU, LDAUTYPE, LDAUL, LDAUPRINT
        <http://cms.mpi.univie.ac.at/vasp/vasp/On_site_Coulomb_interaction_L_S_DA_U.html>`_
  """
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

  def incar_string(self, **kwargs):
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
    return "{0.__class__.__name__}({1!r})".format(self, ["off", "on", "all"][self.value])

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
  def incar_string(self, **kwargs):
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
  def incar_string(self, **kwargs):
    if self.value is None: return None
    return "{0} = {1}".format(self.key.upper(), self.value)
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1}, {2})".format(self, repr(self.key), repr(self.value))

class Choices(SpecialVaspParam):
  """ Vasp parameters with a limited set of choices. 

      Initializes the Choices-type vasp parameters.

      :param key:
          Name of the VASP parameter, e.g. "precfock". It needs not be in
          uppercase. In fact, lower case is preferred for being more pythonic.
      :param choices:
          Dictionary where key is an allowed VASP input for this parameter.
          To each key is associated a list (or set), with allowable forms
          which will translate to the key in the incar. A modified copy of
          this dictionary is owned by the instance being initialized. All
          keys and items should be meaningfully convertible to strings.
      :param default:
          Option from ``choices`` to use as default.

      .. note:: The keys are case-sensitive. The values are not.
  """
  def __init__(self, key, choices, default=None):
    self.key = key
    """ VASP key corresponding to this input. """
    self.choices = {}
    """ Allowable set of choices. """
    for key, items in choices.iteritems():
      self.choices[key] = [u.lower() if hasattr(u, 'lower') else u for u in items]
      self.choices[key].append(key.lower() if hasattr(key, 'lower') else key)
    super(Choices, self).__init__(default)

  @property
  def value(self): return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None; return
    if hasattr(value, 'lower'): value = value.lower()
    for key, items in self.choices.iteritems():
      if value in items: self._value = key; return
    raise ValueError("{0} is not an acceptable choice for {1.key}: {1.choices}.".format(value, self))
  def incar_string(self, **kwargs):
    if self.value is None: return None
    return "{0} = {1}".format(self.key.upper(), self.value)
  def __repr__(self):
    """ Representation of this object. """
    return "{0.__class__.__name__}({1}, {2}, {3})"\
           .format(self, repr(self.key), repr(self.choices), repr(self.value))

class PrecFock(Choices):
  """ Sets up FFT grid in hartree-fock related routines.
      
      Allowable options are:

      - low
      - medium
      - fast
      - normal
      - accurate

      .. note:: The values are not case-sensitive. 
      .. seealso:: `PRECFOCK  <http://cms.mpi.univie.ac.at/vasp/vasp/PRECFOCK_FFT_grid_in_HF_related_routines.html>`_
  """
  def __init__(self, value=None):
    choices = { 'Low': ['low'], 'Medium': ['medium'], 'Fast': ['fast'],
                'Normal': ['normal'], 'Accurate': ['accurate'] }
    super(PrecFock, self).__init__("PRECFOCK", choices, value)
  def __repr__(self):
    return "{0.__class__.__name__}({0.value!r})".format(self)

class Precision(Choices):
  """ Sets accuracy of calculation. 

      - accurate (default)
      - low
      - medium
      - high
      - single

      .. seealso:: `PREC <http://cms.mpi.univie.ac.at/vasp/vasp/PREC_tag.html>`_
  """
  def __init__(self, value = 'accurate'):
    choices = { 'Accurate': ['accurate'], 'Low': ['low'], 'Normal': ['normal'],
                'Medium': ['medium'], 'High': ['high'], 'Single': ['single'] }
    super(Precision, self).__init__('PREC', choices, value)
  def __repr__(self):
    return "{0.__class__.__name__}({0.value!r})".format(self)

class IniWave(Choices):
  """ Specifies how to setup initial wavefunctions.
  
      - 0, jellium
      - 1, random 

      .. seealso:: `INIWAV <http://cms.mpi.univie.ac.at/vasp/guide/node103.html>`_
  """
  def __init__(self, value=None):
    choices = {0: ['jellium'], 1: ['random']}
    super(IniWave, self).__init__('INIWAV', choices, value)
  def __repr__(self):
    return "{0.__class__.__name__}({1!r})".format(self, self.choices[self.value][0])

