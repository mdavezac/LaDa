""" Standard parameter types for use as attributes in Incar """
from ...opt.decorators import broadcast_result
class SpecialVaspParam(object): 
  """ Type checking class. """
  def __init__(self, value): 
    super(SpecialVaspParams, self).__init__()
    self.value = value
  def __repr__(self): return "%s(%s)" % (self.__class__.__name__, self.value)

class NElect(SpecialVaspParam):
  """ Sets number of electrons relative to neutral system.
      
      Gets the number of electrons in the (neutral) system. Then adds value to
      it and computes with the resulting number of electrons.
      >>> nelect = NElect(0) # charge neutral system
      >>> nelect.value = 1   # charge -1 (1 extra electron)
      >>> nelect.value = -1  # charge +1 (1 extra hole)

      @param value: (default:0) number of electrons to add to charge neutral
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
    from boost.mpi import broadcast
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
  def __repr__(self):
    return "%s(%f)" % (self.__class__.__name__, self.value)
          
      
class Algo(SpecialVaspParam): 
  """ Electronic minimization. 
      Can be \"very fast\", \"fast\", or \"normal\" (default). 
  """ 
  def __init__(self, value): super(Algo, self).__init__(value)
  def incar_string(self, *args, **kwargs):
    lower = self.value.lower().rstrip().lstrip()
    lower = lower.replace("_", " ")
    lower = lower.replace("-", " ")
    if lower == "very fast": value = "Very_Fast"
    elif lower == "fast": value = "Fast"
    elif lower == "normal": value = "Normal"
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
  def __init__(self, value): super(Ediff, self).__init__(value)
  def incar_string(self, vasp, *args, **kwargs):
    return "%s = %f " % (self.key, self.value * float(len(vasp._system.atoms)))
  def __repr__(self): return "%s(%e)" % (self.__class__.__name__, repr(self.value))


class Encut(object):
  """ Defines cutoff factor for calculation. 

      There are three ways to set this parameter:
        - if 0 < value <= 3, then the cutoff is value * ENMAX, where ENMAX is
          the maximum recommended cutoff for the species in the system.
        - if value > 3, then prints encut is exactly value.
        - if 0 or None, does not print anything to INCAR
      If 0 or None, uses VASP default.
  """
  def __init__(self, value): super(Ediff, self).__init__(value)

  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    """ Prints ENCUT parameter. """
    from math import fabs
    from ..crystal import specie_list
    assert self.value < -1e-12, ValueError("Wrong value for cutoff.")
    if fabs(self.value) < 1e-12: return "# ENCUT = VASP default"
    if self.value > 3e0 + 1e-12: return "ENCUT = %f" % (self.value)
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
    from ..extract import Extract

    try: return int(self.value[0]), int(self.value[1]), int(self.value[2])
    except IndexError: pass
    except ValueError: pass

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
    from .. import files
    is_root = comm.rank == 0 if comm != None else True
    istart = "0   # start from scratch"
    icharg = "2   # superpositions of atomic densities"
    if self.value == None: istart = "0   # start from scratch"
    elif not self.value.success:
      print "Could not find successful run in directory %s." % (self.value.directory)
      print "Restarting from scratch."
      istart = "0   # start from scratch"
    else:
      ewave = exists( join(self.value.directory, files.WAVECAR) )
      echarge = exists( join(self.value.directory, files.CHGCAR) )
      if ewave:
        path = join(self.value.directory, files.WAVECAR)
        istart = "1  # restart"
        icharg = "0   # from wavefunctions " + path
        if is_root: copy(path, ".")
      elif echarge:
        path = join(self.value.directory, files.CHGCAR)
        istart = "1  # restart"
        icharg = "1   # from charge " + path
        if is_root: copy(path, ".")
      else: 
        istart = "0   # start from scratch"
        icharg = "2   # superpositions of atomic densities"
      if is_root and exists( join(self.value.directory, files.EIGENVALUES) ):
        copy(join(self.value.directory, files.EIGENVALUES), ".") 


    return  "ISTART = %s\nICHARG = %s" % (istart, icharg)

class UParams(SpecialVaspParam): 
  """ Prints U parameters if any found in species settings """
  def __init__(self, *args):
    import re
    
    if verbose == None: value = None
    elif hasattr(verbose, "lower"): 
      verbose = verbose.lower() 
      if verbose == "off": verbose = 0
      elif verbose == "on": verbose = 1
      elif None != re.match(r"\s*occ(upancy)?\s*", verbose): verbose = 1
      elif None != re.match(r"\s*(all|pot(ential)?)\s*", verbose): verbose = 2

    super(UParams, self).__init__(*args)

  @broadcast_result(key=True)
  def incar_string(self, vasp, *args, **kwargs):
    """ Prints LDA+U INCAR parameters. """
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





    
