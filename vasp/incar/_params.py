""" Standard parameter types for use as attributes in Incar """
__docformat__ = "restructuredtext en"
from ...opt.decorators import broadcast_result
class SpecialVaspParam(object): 
  """ Type checking class. """
  def __init__(self, value): 
    super(SpecialVaspParam, self).__init__()
    self.value = value
  def __repr__(self): return "%s(%s)" % (self.__class__.__name__, repr(self.value))

class Magmom(SpecialVaspParam):
  """ Prints magmom to INCAR. 

      There are three types of usage: (i) do nothing, (ii) print a string,
      (iii) use attribute ``magmom`` in the system for which the calculations
      are launched.

      In the second case, the string should be adequately formatted such that
      VASP can understand it for that structure.

      In the third case, the ``magmom`` attribute should be either a (well
      formated) string, or a sequence of floats giving the moment for each
      atom. It is given by setting magmom to the string "magmom", or
      "attribute: somename". This sequence of float must follow the same order
      as the atoms in the system.
  """
  def __init__(self, value = "attribute: magmom"):
    """ Initializes magmom. """
    from re import compile
    super(Magmom, self).__init__(value)

    self._regex = compile("^\s*attribute: (\S+)\s*$")
    """ Regex from which to obtain attribute name. """
    
  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    """ Prints the magmom string if requested. """
    if self.value == None: return None
    if vasp.ispin == 1: return None
    if len(self.value.rstrip().lstrip()) == 0: return None
    if self.value.lower() == "magmom":
      magmom = self._from_attr(vasp, "magmom",  *args, **kwargs)
      return "MAGMOM = {0}".format(magmom) if magmom != None else None
    elif self._regex.match(self.value) != None:
      magmom = self._from_attr(vasp, self._regex.match(self.value).group(1),  *args, **kwargs)
      return "MAGMOM = {0}".format(magmom) if magmom != None else None
    return "MAGMOM = {0}".format(self.value)

  def _from_attr(self, vasp, name, *args, **kwargs):
    """ Creates magmom string from attribute in system. """
    from ...crystal import specie_list
    # checks for existence of the attribute.
    assert hasattr(vasp, "_system"),\
           ValueError("vasp functional does not have a _system attribute.")
    if not hasattr(vasp._system, name): return None
    magmom = getattr(vasp._system, name)
    if magmom == None: return None
    if isinstance(magmom, str): return magmom
    
    # magmom should be a list of moments.
    assert len(magmom) == len(vasp._system.atoms), \
           ValueError("Number of moments and number of atoms does not coincide.")

    result = ""
    for specie in specie_list(vasp._system):
      moments = [u[1] for u in zip(vasp._system.atoms, magmom) if u[0].type == specie]
      tupled = [[1, moments[0]]]
      for m in moments[1:]: 
        if abs(m - tupled[-1][1]) < 1e-12: tupled[-1][0] += 1
        else: tupled.append([1, m])
      for i, m in tupled:
        if i == 1: result += "{0:.2f} ".format(m)
        else:      result += "{0}*{1:.2f} ".format(i, m)
    return result
  
  def __getstate__(self):
    return self.value
  def __setstate__(self, value):
    from re import compile
    self.value = value
    self._regex = compile("^\s*attribute: (\S+)\s*$")

class Npar(SpecialVaspParam):
  """ Parallelization over bands. 

      This parameter is described `here`__.
      It can be set to a particular number:
  
      >>> vasp.npar = 2

      Or it can be deduced automatically. In the latter case, npar is set to the
      largest power of 2 which divides the number of processors:
 
      >>> vasp.npar = "power of two"

      .. __: http://cms.mpi.univie.ac.at/vasp/guide/node138.html
  """

  def __init__(self, value): super(Npar, self).__init__(value)

  def incar_string(self, vasp, *args, **kwargs):
    from re import search
    from math import log
    from ...mpi import Communicator
    if self.value == None: return None
    comm = Communicator(kwargs.get("comm", None))
    if not comm.is_mpi: return None
    if     hasattr(self.value, "lower")\
       and search("power\s+of\s+2", self.value.lower()) != None:
      m = int(log(comm.size)/log(2))
      for i in range(m, -1, -1):
        if comm.size % 2**i == 0: return "NPAR = {0}".format(2**i)
      return None
    else: return "NPAR = {0}".format(self.value)
    



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

  def nelectrons(self, vasp):
    """ Total number of electrons in the system """
    from math import fsum
    # constructs dictionnary of valence charge
    valence = {}
    for key, value in vasp.species.items():
      valence[key] = value.valence
    # sums up charge.
    return fsum( valence[atom.type] for atom in vasp._system.atoms )
    
  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    # gets number of electrons.
    charge_neutral = self.nelectrons(vasp)
    # then prints incar string.
    if self.value == 0:
      return "# NELECT = %s. Charge neutral system" % (charge_neutral)
    elif self.value > 0:
      return "NELECT = %s  # negatively charged system (%i) "\
             % (charge_neutral + self.value, -self.value)
    else: 
      return "NELECT = %s  # positively charged system (+%i) "\
             % (charge_neutral + self.value, -self.value)
          
      
class Algo(SpecialVaspParam): 
  """ Electronic minimization. 
      Can be \"very fast\", \"fast\", or \"normal\" (default). 
  """ 
  def __init__(self, value): super(Algo, self).__init__(value)
  def incar_string(self, *args, **kwargs):
    from .. import is_vasp_5 as _is5
    is_vasp_5 = _is5()
    lower = self.value.lower().rstrip().lstrip()
    lower = lower.replace('_', '')
    lower = lower.replace('-', '')
    if lower == "veryfast": value = "Very_Fast" if not is_vasp_5 else 'VeryFast'
    elif lower in ["fast", 'f']: value = "Fast"
    elif lower in ["normal", 'n']: value = "Normal"
    elif lower in ["damped", 'd']: value = "Damped"
    elif lower == "none" and is_vasp_5: value = "None"
    # below this VASP 5 only options.
    elif not is_vasp_5: raise ValueError, "algo value (%s) is invalid.\n" % (self.value)
    elif lower == "nothing": value = "Nothing"
    elif lower in ["all", 'a']: value = "All"
    elif lower in ["conjugate", 'c']: value = "Conjugate"
    elif lower in ["subrot", 's']: value = "Subrot"
    elif lower in ["eigenval", 'e']: value = "Eigenval"
    elif lower == "chi":  value = "chi"
    elif lower == "gw":   value = "GW"
    elif lower == "gw0":  value = "GW0"
    elif lower == "scgw": value = "scGW"
    elif lower == "scgw0": value = "scGW0"
    else: raise ValueError, "algo value (%s) is invalid.\n" % (self.value)
    return "ALGO = %s" % (value)

class Precision(SpecialVaspParam):
  """ Sets accuracy of calculation. 

      Can be \"accurate\" (default), \"low\", \"medium\", \"high\".
  """
  def __init__(self, value): super(Precision, self).__init__(value)
  def incar_string(self, *args, **kwargs):
    value = self.value.lower()
    if value not in ["accurate", "lower", "medium", "high"]:
      raise ValueError, "PRECISION value (%s) is not allowed." % (self.key)
    return "PRECISION = " + str(self.value)

class Ediff(SpecialVaspParam):
  """ Sets the convergence criteria (per atom) for electronic minimization.

      This tolerance is multiplied by the number of atoms in the system. This
      makes tolerance consistent from one system to the next.
  """
  def __init__(self, value, name="ediff"):
    """ Creates *per atom* tolerance. """
    super(Ediff, self).__init__(value)
    self.name = name
  def incar_string(self, vasp, *args, **kwargs):
    return "{0} = {1} ".format(getattr(self, "name", "ediff").upper(), self.value * float(len(vasp._system.atoms)))
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

  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    """ Prints ENCUT parameter. """
    from math import fabs
    from ...crystal import specie_list
    assert self.value > -1e-12, ValueError("Wrong value for cutoff.")
    if fabs(self.value) < 1e-12: return "# ENCUT = VASP default"
    if self.value > 3e0 + 1e-12:
      return "ENCUT = %f" % (self.value.rescale("eV") if hasattr(self.value, "rescale") else self.value)
    types = specie_list(vasp._system)
    encut = max(vasp.species[type].enmax for type in types)
    return "ENCUT = %f " % (float(encut) * self.value)

class FFTGrid(SpecialVaspParam):
  """ Computes fft grid using VASP. Or if grid is given, computes using that grid. """
  def __init__(self, value): super(FFTGrid, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    return "NGX = %i\nNGY = %i\nNGZ = %i" % self(*args, **kwargs)

  def __call__(self, *args, **kwargs):
    from copy import deepcopy
    from .. import kpoints
    from ...opt.tempdir import Tempdir
    from ...mpi import Communicator
    from ..extract import Extract

    try: return int(self.value[0]), int(self.value[1]), int(self.value[2])
    except IndexError: pass
    except ValueError: pass

    comm = Communicator(kwargs.get('comm', None))
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

class Restart(SpecialVaspParam):
  """ Return from previous run from which to restart.
      
      If None, then starts from scratch.
  """
  def __init__(self, value): super(Restart, self).__init__(value)

  def incar_string(self, *args, **kwargs):
    from os.path import join, exists
    from shutil import copy
    from ...opt import copyfile
    from ...mpi import Communicator
    from .. import files
    comm = Communicator(kwargs.pop("comm", None))
    istart = "0   # start from scratch"
    icharg = "2   # superpositions of atomic densities"
    if self.value == None: istart = "0   # start from scratch"
    elif not self.value.success:
      print "Could not find successful run in directory %s." % (self.value.directory)
      print "Restarting from scratch."
      istart = "0   # start from scratch"
    else:
      comm.barrier()
      ewave = exists( join(self.value.directory, files.WAVECAR) )
      echarge = exists( join(self.value.directory, files.CHGCAR) )
      if ewave:
        path = join(self.value.directory, files.WAVECAR)
        istart = "1  # restart"
        icharg = "0   # from wavefunctions " + path
        if comm.is_root: copy(path, ".")
      elif echarge:
        path = join(self.value.directory, files.CHGCAR)
        istart = "1  # restart"
        icharg = "1   # from charge " + path
        if comm.is_root: copy(path, ".")
      else: 
        istart = "0   # start from scratch"
        icharg = "2   # superpositions of atomic densities"
      if comm.is_root:
        copyfile(join(self.value.directory, files.EIGENVALUES), nothrow='same exists') 
        copyfile(join(self.value.directory, files.CONTCAR), files.POSCAR,
                 nothrow='same exists', symlink=getattr(args[0], 'symlink', False)) 
      comm.barrier()
    return  "ISTART = %s\nICHARG = %s" % (istart, icharg)

class UParams(SpecialVaspParam): 
  """ Prints U parameters if any found in species settings """
  def __init__(self, value):
    import re
    
    if value == None: value = None
    elif hasattr(value, "lower"): 
      value = value.lower() 
      if value == "off": value = 0
      elif value == "on": value = 1
      elif None != re.match(r"\s*occ(upancy)?\s*", value): value = 1
      elif None != re.match(r"\s*(all|pot(ential)?)\s*", value): value = 2

    super(UParams, self).__init__(value)

  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    """ Prints LDA+U INCAR parameters. """
    from ...crystal import specie_list
    types = specie_list(vasp._system)
    # existence and sanity check
    has_U, which_type = False, None 
    for type in types:
      specie = vasp.species[type]
      if len(specie.U) == 0: continue
      if len(specie.U) > 4: 
        raise AssertionError, "More than 4 channels for U/NLEP parameters"
      has_U = True
      # checks consistency.
      which_type = specie.U[0]["type"]
      for l in specie.U[1:]: 
        assert which_type == l["type"], \
               AssertionError("LDA+U/NLEP types are not consistent across species.")
    if not has_U: return "# no LDA+U/NLEP parameters ";

    # Prints LDA + U parameters
    result = "LDAU = .TRUE.\nLDAUPRINT = %i\nLDAUTYPE = %i\n" % (self.value, which_type)

    for i in range( max(len(vasp.species[type].U) for type in types) ):
      line = "LDUL%i=" % (i+1), "LDUU%i=" % (i+1), "LDUJ%i=" % (i+1), "LDUO%i=" % (i+1)
      for type in types:
        specie = vasp.species[type]
        a = -1, 0e0, 0e0, 1
        if len(specie.U) <= i: pass
        elif specie.U[i]["func"] == "U":    
          a = [specie.U[i]["l"], specie.U[i]["U"], specie.U[i]["J"], 1]
        elif specie.U[i]["func"] == "nlep": 
          a = [specie.U[i]["l"], specie.U[i]["U"], 0e0, 2]
        elif specie.U[i]["func"] == "enlep":
          a = [specie.U[i]["l"], specie.U[i]["U0"], specie.U[i]["U1"], 3]
        else: raise RuntimeError, "Debug Error."
        if hasattr(a[1], "rescale"): a[1] = a[1].rescale("eV")
        if hasattr(a[2], "rescale"): a[2] = a[2].rescale("eV")
        line = "%s %i"      % (line[0], a[0]),\
               "%s %18.10e" % (line[1], a[1]),\
               "%s %18.10e" % (line[2], a[2]),\
               "%s %i"      % (line[3], a[3])
      result += "\n%s\n%s\n%s\n%s\n" % (line)
    return result

class IniWave(SpecialVaspParam):
  def __init__(self, value): super(IniWave, self).__init__(value)
  @broadcast_result(key=True)
  def incar_string(self, *args, **kwargs):
    """ Returns VASP incar string. """
    if self.value == "1" or self.value == "random": result = 1
    elif self.value == "0" or self.value == "jellium": result = 0
    else: raise ValueError("iniwave cannot be set to " + self.value + ".")
    return "INIWAVE = %i\n"  % (value)


