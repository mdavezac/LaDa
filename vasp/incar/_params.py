""" Standard parameter types for use as attributes in Incar """
from ...opt.decorators import broadcast_result as _bcastresult
class Standard(object):
  """
      A standard parameter in the form of key/value pair. 
      
      In practice, the key should be an INCAR tag, and the value something
      which makes sense for that key.
      C{print self.incar_string(vasp, comm=comm)} outputs C{print
      str(self.key), "=", str(self.value)}.
  """
  def __init__(self, key, value, validity=None):
    super(Standard, self).__init__()

    self.key = key
    """ Any string """
    self._value = value
    """ Any value """
    self.validity = validity
    """ Functor returning true if the value is valid. """

  # getter setter functions for value.
  def _getvalue(self): return self._value
  def _setvalue(self, value): 
    if self.validity != None:
      if not self.validity(value):
        raise ValueError, "Value %s = %s is invalid.\n%s\n"  % (self.key, value, self.__doc__)
    self._value = value

  value = property(_getvalue, _setvalue)
  """ Actual property 
  
      It takes a meaningful value when set (eg a string), and
      returns its  "VASP" equivalent when gotten.
  """

  @_bcastresult
  def incar_string(self, vasp):
    """ Prints the key/value as key = value """
    return "%s = %s" % (self.key, str(self.value))

  def __getstate__(self):
    from marshal import dumps
    d = self.__dict__.copy()
    del d["validity"]
    if self.validity == None: return d, None
    return d, dumps(self.validity.func_code)
  def __setstate__(self, arg):
    from marshal import loads
    self.__dict__.update(arg[0])
    self.validity = None
    if arg[1] != None: 
      self.validity = lambda x: x
      self.validity.func_code  = loads(arg[1])


class NoPrintStandard(Standard):
  """ Does not print if value is the default given on initialization. """
  def __init__(self, *args, **kwargs):
    super(NoPrintStandard, self).__init__(*args, **kwargs)
    self.default = self.value

  @_bcastresult
  def incar_string(self, vasp):
    if self.default == self.value: 
      return "# %s = VASP default." % (self.key)
    return super(NoPrintStandard, self).incar_string(vasp)
  
  def __setstate__(self, arg): 
    super(NoPrintStandard, self).__setstate__(arg)
    self.default = self.value

class Iniwave(Standard):
  """ Iniwave variable 

      self.value can be either "random" or "jellium".
  """
  def __init__(self, value = "random"):
    validity = lambda x: x=="random" or x=="jellium"
    super(Iniwave, self).__init__("INIWAVE", value, validity)

  def _getvalue(self): 
    if self._value == "random": return 1
    elif self._value == "jellium": return 0
    else: raise RuntimeError, "Internal bug."

  value = property(_getvalue, Standard._setvalue)
  """ Iniwave property
  
      It takes a meaningful value when set (C{self.value = random})
        - random: wavefunctions are initalized with random number. Safest. 
        - jellium: wavefunctions are initalized a jelliums. Not so safe.
      returns its  "VASP" equivalent when gotten (C{print self.value} outputs 0 or 1).
  """


class Algo(Standard): 
  """ Electronic minimization. 
      Can be \"very fast\", \"fast\", or \"normal\" (default). 
  """ 
  def __init__(self):
    super(Algo, self).__init__("ALGO", "normal") 
  def _getvalue(self):
    if self._value not in ["Very_Fast", "Fast", "Normal"]:
      self.value = self._value
    return self._value
  def _setvalue(self, value):
    lower = value.lower().rstrip().lstrip()
    lower = lower.replace("_", " ")
    lower = lower.replace("-", " ")
    if lower == "very fast": self._value = "Very_Fast"
    elif lower == "fast": self._value = "Fast"
    elif lower == "normal": self._value = "Normal"
    else: raise ValueError, "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)
  value = property(_getvalue, _setvalue)

class Precision(Standard):
  """ Sets accuracy of calculation. 
      Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
  """
  def __init__(self): 
    super(Precision, self).__init__("PREC", "accurate") 
  def _getvalue(self):
    if self._value not in ["accurate", "lower", "medium", "high"]:
      self.value = self._value
    return self._value
  def _setvalue(self, value):
    self._value = value.lower()
    if self._value not in ["accurate", "lower", "medium", "high"]:
      raise ValueError, "%s = %s is invalid.\n%s\n" % (self.key, value, self.__doc__)
  value = property(_getvalue, _setvalue)

class Ediff(Standard):
  """ Sets the convergence criteria (per atom) for electronic minimization.

      This tolerance is divided by the number of atoms in the system. This
      makes tolerance consistent from one system to the next.
  """
  def __init__(self, value = 1e-4):
    super(Ediff, self).__init__("EDIFF", value, validity = lambda x: x > 0e0)

  @_bcastresult
  def incar_string(self, vasp):
    return "%s = %f " % (self.key, self.value * float(len(vasp._system.atoms)))


class Encut(object):
  """ Defines energy cutoff for calculation. 
      The base cutoff is obtained from the species of the owner specified in
      __init__. The actual cutoff is that times the safety (also specified in __init__).
  """
  class Safety(object):
    def __init__(self, value): self.value = value
    def __set__(self, owner, value):
      assert value > 0, "Ecut safety must be larger than 0.\n"
      self._value = value
    def __get__(self, owner, parenttype=None): return self._value

  safety = Safety(1.25)
  """ A safety upon the largest ENMAX """
  key = "ENCUT"
  """ INCAR key """

  def __init__(self, safety = 1.25):
    super(Encut, self).__init__()
    self.safety = safety

  def _getvalue(self):
    raise ValueError, "Value of encut cannot be obtained directly. use enmax method instead.\n"
  def _setvalue(self, value):
    raise ValueError, "Value of encut cannot be set manually. Change safety instead.\n"
  x = property(_getvalue, _setvalue)

  @_bcastresult
  def incar_string(self, vasp):
    from math import fabs
    if fabs(self.safety - 1e0) < 1e-12: return "# ENCUT = VASP default"
    return "%s = %f " % (self.key, float(self.enmax(vasp.species)) * self.safety)

  @staticmethod
  def enmax(species):
    """ Retrieves max ENMAX from list of species. """
    from math import ceil
    return ceil( max(s.enmax for s in species) )

class Smearing(object):
  """ Value of the smearing used in the calculation. 
      It can be specified as a string: "type x", where type is any of fermi,
      gaussian, mp, tetra, metal or insulator, and x is the energy scale in eV.
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
  def __init__(self, string = "insulator 0.2"):
    super(Smearing, self).__init__()
    self.value = string

  def _getvalue(self): return (self.type, self.smearing)
  def _setvalue(self, value):
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
      if second != None: self.smearing = float(second)
    elif first == "gaussian" or first == "0":
      self.type == 0
      if second != None: self.smearing = float(second)
    elif first == "metal":
      self.type = 1
      if second != None: self.smearing = float(second)
    elif first == "mp" or first == "metal":
      if third != None:
        self.type = int(second)
        self.smearing = float(third)
      elif second != None:
        self.type = int(second)
      else: self.type = 1
      assert self.type >= 1, "Mehtfessel-Paxton order must be at least 1."
    elif first == "bloechl" or first == "-5" or first == "insulator":
      self.type = -5
      if second != None: self.smearing = float(second)
    elif first == "tetra" or first == "-4":
      self.type = -4
      if second != None: self.smearing = float(second)
    else: 
      try: self._value = int(first)
      except: raise ValueError, "Unknown smearing value %s.\nsmearing: %s\n" % (value, self.__doc__)
      assert self._value >= 1, "Unknown smearing value %s.\nsmearing: %s\n" % (value, self.__doc__)
  value = property(_getvalue, _setvalue)
  """
      Type of smearing used.
      If "gotten", this decorator returns a tuple containing the type and the energy scale.
      It can be set using a similar tuple, or a string input with keywords. 
  """
  @_bcastresult
  def incar_string(self, vasp):
    return "ISMEAR  = %s\nSIGMA   = %f\n" % self.value

  
class Isym(Standard):
  """ Type of symmetry used in the calculation.
      Can be "off" or a float corresponding to the tolerance used to determine
      symmetry operation.  By default, it is 1e-5.
  """
  def __init__(self, value = 1e-5):
    super(Isym, self).__init__("SYMPREC", value, None)
    self._value = value

  def _getvalue(self):
    if self._value <= 0e0: return "off"
    return self._value
  def _setvalue(self, value):
    strval = str(value)
    if strval == "off": self._value = -1e0
    else: self._value = float(value) 
  value = property(_getvalue, _setvalue)

  @_bcastresult
  def incar_string(self, vasp):
    if self._value <= 0e0: return "ISYM = 0\n" 
    return super(Isym,self).incar_string(vasp)

  def __eq__(self, other):
    if isinstance(other, Isym):
      return self._value == self._value
    return self == Isym(other)
  
class FFTGrid(object):
  """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """
  def __init__(self, grid = None): 
    super(FFTGrid, self).__init__()
    self._value = grid

  def _getvalue(self):
    if self._value == None:
      raise "Value is computed from VASP. Call self.__call__(vasp) instead.\n"
    return self._value

  def _setvalue(self, value):
    self._value = value
    if self._value != None: 
      assert len(self._value) == 3, "FFT grid should be none or a 3-tuple.\n"
      assert self._value[0] > 0, "FFT grid components should be positive.\n"
      assert self._value[1] > 0, "FFT grid components should be positive.\n"
      assert self._value[2] > 0, "FFT grid components should be positive.\n"
  value = property(_getvalue, _setvalue)

  def incar_string(self, vasp, comm=None):
    return "NGX = %i\nNGY = %i\nNGZ = %i" % self(vasp, comm)

  def __call__(self, vasp, comm):
    from copy import deepcopy
    from os import getcwd
    from .. import kpoints
    from ...opt.tempdir import Tempdir
    from ..extract import Extract

    if self._value != None: return self._value

    vasp = deepcopy(vasp)
    with Tempdir(workdir=vasp._tempdir, keep=True, comm=comm) as vasp._tempdir:
      vasp.kpoints = kpoints.Gamma()
      vasp.relaxation.value = None
      # may need to do some checking here... will run into infinit recursion
      # if attribute is not fftgrid. Can we access vasp.__dict__ and remove
      # reference to self?
      del vasp.fftgrid # simply remove attribute to get VASP default
      # makes sure we do nothing during this run
      vasp.nelmdl = Standard("NELMDL",  0)
      vasp.nelm   = Standard("NELM",  0)
      vasp.nelmin = Standard("NELMIN",  0)

      # Now runs vasp. OUTCAR should be in temp indir
      vasp._prerun(comm)
      vasp._run(comm)
      # no need for postrun

      # finally extracts from OUTCAR.
      return Extract(directory = vasp._tempdir, comm = comm).fft



      

class Restart(object):
  """
      Directory where to restart, or None.
      
      If None, then starts from scratch.
      If this directory contains WAVECAR or CHGCAR, then restarts from
      wavefunctions or charge. If this directory does not exist, issue an
      error.
  """
  def __init__(self, directory = None):  
    super(Restart, self).__init__()
    self._value = directory

  def __get__(self, owner, ownertype=None): return self
  def __set__(self, owner, value): self._setvalue(self, value)
  def _getvalue(self): return self._value
  def _setvalue(self, object): self._value = object
  value = property(_getvalue, _setvalue)

  @_bcastresult
  def incar_string(self, vasp):
    istart = "0   # start from scratch"
    icharg = "2   # superpositions of atomic densities"
    if self.value == None or len(self.value) == 0:
      istart = "0   # start from scratch"
    elif not os.path.exists( self.value ):
      raise RuntimeError, "Could not find restart directory " + self.value + "\n";
    else:
      ewave = os.path.exists( os.path.join( self.value, 'WAVECAR' ) )
      echarge = os.path.exists( os.path.join( self.value, 'CHGCAR' ) )
      if ewave:
        istart = "1  # restart"
        icharg = "0   # from wavefunctions " + os.path.join( self.value, 'WAVECAR' )
      elif echarge:
        istart = "1  # restart"
        icharg = "1   # from charge " + os.path.join( self.value, 'CHGCAR' )
      else: 
        istart = "0   # start from scratch"
        icharg = "2   # superpositions of atomic densities"

    return  "ISTART = %s\nICHARG = %s" % (istart, icharg)

class Relaxation(object):
  """ Sets ISIF in incar depending on type relaxation required. 

      @param value: Specifies the type of relaxation. Defaults to None. 
        - if set to None or empty string, then no relaxation.
        - if ionic is in string, then includes ionic relaxation.
        - if cellshape is in string, then includes cell relaxation.
        - if volume is in string, then includes volume relaxation.
        Makes sure that the parameters make sense together. 
        Can also be set using an integer between 0 and 7. See VASP manual. 
      @param nsw: the number of ionic steps to perform. Defaults to 40.
      @param ibrion: the type of ionic relaxation algorithm. Defaults to 2.
      @param potim: ionic-time step during relaxation. Defaults to 0.5.
  """
  key = "ISIF"
  """ INCAR key """

  def __init__(self, value = None, nsw=40, ibrion=2, potim=0.5):
    super(Relaxation, self).__init__()
    self._value = value
    self.nsw = 40
    self.ibrion = 2
    self.potim = 0.5

  def _getvalue(self): return self._value
  def _setvalue(self, object):
    import re

    if object == None:
      self._value = 0
      return 

    try:
      self._value = int(object) 
      assert self._value >= 0
      return 
    except ValueError: pass

    if object.lower() == "static": 
      self._value = 0
      return 

    ionic =  re.search( "ion(ic|s)?", object.lower() ) != None
    cellshape = re.search( "cell(\s+|-|_)?(?:shape)?", object.lower() ) != None
    volume = re.search( "volume", object.lower() ) != None
    if (not ionic) and (not cellshape) and (not volume): self._value = 0
    elif ionic and (not cellshape) and (not volume):     self._value = 2
    elif ionic and cellshape and volume:                 self._value = 3
    elif ionic and cellshape and (not volume):           self._value = 4
    elif(not ionic) and  cellshape and (not volume):     self._value = 5
    elif(not ionic) and  cellshape and volume:           self._value = 6
    elif(not ionic) and (not cellshape) and volume:      self._value = 7
    elif ionic and (not cellshape) and volume:
      raise RuntimeError, "VASP does not allow relaxation of atomic position"\
                          "and volume at constant cell-shape.\n"
  value = property(_getvalue, _setvalue)

  @_bcastresult
  def incar_string(self, vasp):
    isif = self.value
    result = "NSW    = %3i   # number of ionic steps.\n"\
             "IBRION =   %1i # ionic-relaxation minimization method.\n"\
             "POTIM  = %7.4f # ionic-relaxation step \n" % (self.nsw, self.ibrion, self.potim)
    if isif == 0: return "ISIF = 0     # static calculation " 
    elif isif == 1:       
      return result + "ISIF   = 2     # relaxing atomic positions. Only tr(strain) is correct."
    elif isif == 2:           
      return result + "ISIF   = 2     # relaxing atomic positions."
    elif isif == 3:           
      return result + "ISIF   = 3     # relaxing all structural degrees of freedom."
    elif isif == 4:           
      return result + "ISIF   = 4     # relaxing atom. pos. and cell-shape at constant V."
    elif isif == 5:           
      return result + "ISIF   = 5     # relaxing cell-shape at constant atom.-pos. and V."
    elif isif == 6:           
      return result + "ISIF   = 6     # relaxing V. and cell-shape at constant atom.-pos."
    elif isif == 7:           
      return result + "ISIF   = 7     # relaxing V. only."
    raise RuntimeError("Internal bug.")

class NElect(object):
  """ Sets number of electrons relative to neutral system.
      
      Gets the number of electrons in the (neutral) system. Then adds value to
      it and computes with the resulting number of electrons.
      >>> nelect = NElect(0) # charge neutral system
      >>> nelect.value = 1   # charge -1 (1 extra electron)
      >>> nelect.value -= 2  # charge +1 (1 less electron)

      @param value: (default:0) number of electrons to add to charge neutral
                    system.
  """
  key = "NELECT"
  """ INCAR Key """

  def __init__(self, value = 0):
    super(NElect, self).__init__()
    self.value = value

  def nelectrons(self, vasp):
    """ Total number of electrons in the system """
    from math import fsum
    # constructs dictionnary of valence charge
    valence = {}
    for s in vasp.species:
      valence[s.symbol] = s.valence
    # sums up charge.
    return fsum( valence[atom.type] for atom in vasp._system.atoms )
    
  @_bcastresult
  def incar_string(self, vasp):
    from boost.mpi import broadcast
    # gets number of electrons.
    charge_neutral = self.nelectrons(vasp)
    # then prints incar string.
    if self.value == 0:
      return "# %s = %s. Charge neutral system" % (self.key, charge_neutral)
    elif self.value > 0:
      return "%s = %s  # negatively charged system (%i) "\
             % (self.key, charge_neutral + self.value, -self.value)
    else: 
      return "%s = %s  # positively charged system (+%i) "\
             % (self.key, charge_neutral + self.value, -self.value)
          
class NBands(NElect):
  """ Sets number of bands to compute relative to number of electrons. 
      
      Three types of input are accepted:
  """
  default = -1
  """ Use vasp default """
  nunit = 0
  """ Adds n bands, with n=self.value. """
  nions = 1
  """ Adds n bands, with n=self.value * \"number of atoms in system\". """
  nelec = 2
  """ Adds n bands, with n=self.value * \"number of electrons in system\". """
  def __init__(self, value = 0, unit = None):
    super(NBands, self).__init__(value)
    if unit == None: self.unit = default
    else:            self.unit = unit

  key = "NBANDS"
  """ INCAR key """

  @_bcastresult
  def incar_string(self, vasp):
    from boost.mpi import broadcast
    # gets number of electrons.
    ne = self.nelectrons(vasp)
    # returns adequate string.
    if self.unit == NBands.default:
      return "# %s = VASP default " % (self.key)
    elif self.unit == NBands.nunit:
      return "%s = %s" % (self.key, self.value + float(ne))
    elif self.unit == NBands.nions:
      return "%s = %s" % (self.key, int(ne + ceil(self.value * float(len(vasp._system.atoms)))) )
    elif self.unit == NBands.nelec:
      return "%s = %s" % (self.key, int(ne + ceil(self.value * float(ne))) )
    else: raise "Unknown method in NBands"
      
class UParams(object): 
  """ Prints U parameters if any found in species settings """
  def __init__(self, verbose=None, **kwargs):
    import re
    
    if verbose == None: self.value = 0
    elif hasattr(verbose, "lower"): 
      verbose = verbose.lower() 
      if verbose == "off": verbose = 0
      elif verbose == "on": verbose = 1
      elif None != re.match(r"\s*occ(upancy)?\s*", verbose): verbose = 1
      elif None != re.match(r"\s*(all|pot(ential)?)\s*", verbose): verbose = 2

    self.value = int(verbose)
    super(UParams, self).__init__(**kwargs)

  @_bcastresult
  def incar_string(self, vasp):
    # existence and sanity check
    has_U, which_type = False, None 
    for specie in vasp.species:
      if len(specie.U) == 0: continue
      if len(specie.U) > 4: 
        raise AssertionError, "More than 4 channels for U/NLEP parameters"
      has_U = True
      for l in specie.U: 
        if which_type == None:
          which_type = l["type"]
        elif which_type != l["type"]:
          raise AssertionError, "LDA+U/NLEP types are not consistent across species."
    if not has_U: return "# no LDA+U/NLEP parameters ";

    result = "LDAU = .TRUE.\nLDAUPRINT = %i\nLDAUTYPE = %i\n" % (self.value, which_type)

    for i in range( max(len(specie.U) for specie in vasp.species) ):
      line = "LDUL%i=" % (i+1), "LDUU%i=" % (i+1), "LDUJ%i=" % (i+1), "LDUO%i=" % (i+1)
      for specie in vasp.species:
        a = -1, 0e0, 0e0, 1
        if len(specie.U) <= i: pass
        elif specie.U[i]["func"] == "U":    
          a = specie.U[i]["l"], specie.U[i]["U"], specie.U[i]["J"], 1
        elif specie.U[i]["func"] == "nlep": 
          a = specie.U[i]["l"], specie.U[i]["U"], 0e0, 2
        elif specie.U[i]["func"] == "enlep":
          a = specie.U[i]["l"], specie.U[i]["U0"], specie.U[i]["U1"], 3
        else: raise RuntimeError, "Debug Error."
        line = "%s %i"      % (line[0], a[0]),\
               "%s %18.10e" % (line[1], a[1]),\
               "%s %18.10e" % (line[2], a[2]),\
               "%s %i"      % (line[3], a[3])
      result += "\n%s\n%s\n%s\n%s\n" % (line)
    return result





    
