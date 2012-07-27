from ..functools.keyword import BoolKeyword as BaseBoolKeyword, ValueKeyword,  \
                                TypedKeyword, AliasKeyword
class BoolKeyword(BaseBoolKeyword):
  """ Boolean keyword.

      If True, the keyword is present.
      If False, it is not.
  """
  def __init__(self, keyword=None, value=False):
    """ Initializes FullOptG keyword. """
    super(BoolKeyword, self).__init__(keyword=keyword, value=value)
  def input_map(self, **kwargs):
    """ Map keyword, value """
    if self.value is None: return None
    if getattr(self, 'keyword', None) is None: return None
    return { self.keyword: '.TRUE.' if self.value else '.FALSE.' }

class Magmom(ValueKeyword):
  """ Sets the initial magnetic moments on each atom.

      There are three types of usage: 

      - if None or False, does nothing
      - if calculations are not spin-polarized, does nothing.
      - if a string, uses that as for the MAGMOM_ keyword
      - if True and at least one atom in the structure has a non-zero
        ``magmom`` attribute, then creates the relevant moment input for VASP_

      If the calculation is **not** spin-polarized, then the magnetic moment
      tag is not set.

      .. note:: Please set by hand for non-collinear calculations

      .. seealso:: MAGMOM_

      .. _MAGMOM: http://cms.mpi.univie.ac.at/wiki/index.php/MAGMOM
  """
  keyword = 'MAGMOM'
  'VASP keyword'
  def __init__(self, value=True):
    super(Magmom, self).__init__(value)

  def value(self):
    """ MAGMOM value, or whether to compute it. 

        - if None or False, does nothing
        - if calculations are not spin-polarized, does nothing.
        - if a string, uses that as for the MAGMOM_ keyword
        - if True and at least one atom in the structure has a non-zero
          ``magmom`` attribute, then creates the relevant moment input for
          VASP_
    """
    return self._value
  @value.setter
  def value(self, value):
    if value is None: self._value = None
    elif value is True or value is False: self._value = value
    elif not isinstance(value, str):
      raise ValueError( 'Unknown value for magmom {0}. '                       \
                        'Should be True, False, None, or an adequate string'   \
                        .format(value) )
    else: self._value = value

  def output_map(self, **kwargs):
    """ MAGMOM input for VASP. """
    from ...crystal import specieset
    if self.value is None or self.value == False: return None
    if kwargs['vasp'].ispin == 1: return None
    if isinstance(self.value, str): return {self.keyword: str(self.value)}
    
    structure = kwargs['structure']
    if all(not hasattr(u, 'magmom') for u in structure): return None
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
    return { self.keyword, result.rstrip() }

class System(ValueKeyword):
  """ System title to use for calculation.

      - If None and ... 
         - if the structure has a ``name`` attribute, uses that as the
           calculations title
         - else does not use SYSTEM_ tag
      - If something else which is convertable to a string,  and ...
         - if the structure has a ``name`` attribute, uses 'string: name' as
           the title
         - otherwise, uses the string

      .. seealso:: SYSTEM_

      .. _SYSTEM: http://cms.mpi.univie.ac.at/vasp/guide/node94.html>
  """
  keyword = 'SYSTEM'
  """ VASP keyword """
  def __init__(self, value): super(System, self).__init__(value=value)

  def output_map(self, **kwargs):
    """ Tries to return sensible title. 

        Never throws.
    """
    try: 
      if self.value is None:
        if len(getattr(kwargs['structure'], 'name', '')) == 0: return None
        return { self.keyword: str(kwargs['structure'].name) }
      elif len(getattr(kwargs['structure'], 'name', '')) == 0:
        try: return { self.keyword: str(self.value) }
        except: return None
      try:
        return { self.keyword:
                 '{0}: {1}'.format(str(self.value), kwargs['structure'].name) }
      except: return { self.keyword: str(kwargs['structure'].name) }
    except: return None

class Npar(ValueKeyword):
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
  keyword = 'NPAR'
  """ VASP keyword """
  def __init__(self, value): super(Npar, self).__init__(value=value)

  def output_map(self, **kwargs):
    from math import log, sqrt
    if self.value is None: return None
    if not isinstance(self.value, str): 
      if self.value < 1: return None
      return {self.keyword: str(self.value)}

    comm = kwargs.get('comm', None)
    if comm is None: return None
    n = comm['n']
    if n == 1: return None
    if self.value == "power of two":
      m = int(log(n)/log(2))
      for i in range(m, 0, -1):
        if n % 2**i == 0: return {self.keyword: str(i)}
      return None
    if self.value == "sqrt":
      return {self.keyword: str(int(sqrt(n)+0.001))}
    raise ValueError("Unknown request npar = {0}".format(self.value))

class ExtraElectron(TypedKeyword):
  """ Sets number of electrons relative to neutral system.
      
      Gets the number of electrons in the (neutral) system. Then adds value to
      it and computes with the resulting number of electrons.

      >>> vasp.extraelectron =  0  # charge neutral system
      >>> vasp.extraelectron =  1  # charge -1 (1 extra electron)
      >>> vasp.extraelectron = -1  # charge +1 (1 extra hole)

      Disables :py:attr:`lada.vasp.functional.Functional.nelect` if set to
      something other than None.

      .. seealso:: `NELECT <http://cms.mpi.univie.ac.at/wiki/index.php/NELECT>`_
  """
  type = float
  """ Type of this input. """
  keyword = 'NELECT'
  """ VASP keyword. """
  def __init__(self, value=0): super(ExtraElectron, self).__init__(value=value)

  def __set__(self, instance, value):
    if value is not None: instance.nelect = None
    return super(ExtraElectron, self).__set__(instance, value)

  def nelectrons(self, vasp, structure):
    """ Total number of electrons in the system """
    from math import fsum
    # constructs dictionnary of valence charge
    valence = {}
    for key, value in vasp.species.items():
      valence[key] = value.valence
    # sums up charge.
    return fsum( valence[atom.type] for atom in structure )
    
  def output_map(self, **kwargs):
    # gets number of electrons.
    charge_neutral = self.nelectrons(kwargs['vasp'], kwargs['structure'])
    # then prints incar string.
    if self.value == 0: return None
    return {self.keyword: str(charge_neutral + self.value)}

class NElect(TypedKeyword):
  """ Sets the absolute number of electrons.
      
      Disables :py:attr:`lada.vasp.functional.Functional.extraelectron` if set to
      something other than None.

      .. seealso:: `NELECT <http://cms.mpi.univie.ac.at/wiki/index.php/NELECT>`_
  """
  type = float
  """ Type of this input. """
  keyword = 'NELECT'
  """ VASP keyword. """
  def __init__(self, value=0): super(ExtraElectron, self).__init__(value=value)

  def __set__(self, instance, value):
    if value is not None: instance.extraelectron = None
    return super(NElect, self).__set__(instance, value)

class Algo(ValueKeyword): 
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
  keyword = 'ALGO'
  """ VASP keyword. """
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
    if is_vasp_4                                                               \
       and ( lower[0] in ['c', 's', 'e']                                       \
             or lower in [ "nothing", "subrot", "exact",                       \
                           "gw", "gw0", "chi", "scgw",                         \
                           "scgw0"] ): 
      raise ValueError("algo value ({0}) is not valid with VASP 4.6.".format(value))

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
      raise ValueError("algo value ({0!r}) is invalid.".format(value))
    self._value = value

class Ediff(TypedKeyword):
  """ Sets the convergence criteria (per atom) for electronic minimization.

      - value > 0e0: the tolerance is multiplied by the number of atoms in the
        system. This makes tolerance consistent from one system to the next.
      - value < 0e0: tolerance is given as absolute value, without multiplying
        by size of system.

      .. seealso:: `EDIFF <http://cms.mpi.univie.ac.at/vasp/guide/node105.html>`_
  """
  type = float
  """ Type of the value """
  keyword = 'EDIFF'
  """ VASP keyword """
  def __init__(self, value=1e0):
    """ Creates *per atom* tolerance. """
    super(Ediff, self).__init__(value=value)
  def output_map(self, **kwargs):
    if self.value is None: return 
    if self.value < 0: return {self.keyword: str(-self.value)}
    return { self.keyword: str(self.value * float(len(kwargs["structure"]))) }

class Ediffg(Ediff):
  """ Sets the convergence criteria (per atom) for ionic minimization.

      - value > 0e0: the tolerance is multiplied by the number of atoms in the
        system. This makes tolerance consistent from one system to the next.
      - value < 0e0: tolerance is given as is (negative), and applies to forces.

      .. seealso:: `EDIFFG <http://cms.mpi.univie.ac.at/vasp/guide/node107.html>`_
  """
  keyword = 'EDIFFG'
  def __init__(self, value=None):
    """ Creates *per atom* tolerance. """
    super(Ediffg, self).__init__(value)
  def output_map(self, **kwargs):
    if self.value is None: return 
    if self.value < 0: return {self.keyword: str(self.value)}
    return { self.keyword: str(self.value * float(len(kwargs["structure"]))) }

class Encut(ValueKeyword):
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
  keyword = "ENCUT"
  """ Corresponding VASP key. """
  def __init__(self, value=None): super(Encut, self).__init__(value=value)

  def output_map(self, **kwargs):
    from quantities import eV
    from ...crystal import specieset
    value = self._value
    if hasattr(self.value, 'units'):
      value = self.value.rescale(eV).magnitude
      return {self.keyword: str(value)} if value > 1e-12 else None
    if value is None:   return None
    elif value < 1e-12: return None
    elif value >= 1e-12 and value <= 3.0:
      types = specieset(kwargs["structure"])
      encut = max(kwargs["vasp"].species[type].enmax for type in types)
      if hasattr(encut, 'rescale'): encut = float(encut.rescale(eV))
      return {self.keyword: str(encut * value)}
    return {self.keyword: str(value)}

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
  keyword = 'ENCUTGW'

class ICharge(ValueKeyword):
  """ Charge from which to start. 

      It is best to keep this attribute set to -1, in which case, LaDa takes
      care of copying the relevant files.

        - -1: Automatically determined by LaDA. Depends on the value of restart_
              and the existence of the relevant files. Also takes care of non-scf
              bit.
  
        - 0: Tries to restart from wavefunctions. Uses the latest WAVECAR file
             between the one currently in the output directory and the one in
             the restart directory (if speciefied). Sets nonscf_ to False.
  
             .. note:: CHGCAR is also copied, just in case.
  
        - 1: Tries to restart from wavefunctions. Uses the latest WAVECAR file
             between the one currently in the output directory and the one in
             the restart directory (if speciefied). Sets nonscf_ to False.
  
        - 2: Superimposition of atomic charge densities. Sets nonscf_ to False.
  
        - 4: Reads potential from POT file (VASP-5.1 only). The POT file is
             deduced the same way as for CHGAR and WAVECAR above.  Sets nonscf_
             to False.
  
        - 10, 11, 12: Same as 0, 1, 2 above, but also sets nonscf_ to True. This
             is a shortcut. The value is actually kept to 0, 1, or 2:
  
             >>> vasp.icharg = 10
             >>> vasp.nonscf, vasp.icharg
             (True, 0)

      .. note::
      
         Files are copied right before the calculation takes place, not before.

      .. seealso:: ICHARG_

      .. _ICHARG: http://cms.mpi.univie.ac.at/wiki/index.php/ICHARG
      .. _restart: :py:attr:`~lada.vasp.functional.Functional.restart`
      .. _nonscf: :py:attr:`~lada.vasp.functional.Functional.nonscf`
  """ 
  keyword = 'ICHARG'
  """ VASP keyword """
  def __init__(self, value=-1): 
    super(ICharge, self).__init__(value)
  def __set__(self, instance, value):
    """ Sets internal value. 

        Makes sure that the input value is allowed, and that the nonscf_
        attribute is set properly.

      .. _nonscf: :py:attr:`~lada.vasp.functional.Functional.nonscf`
    """
    from ..error import ValueError
    if value is None: 
      self.value = None
      return
    value = int(value)
    if value not in [-1, 0, 1, 2, 4, 10, 11, 12]: 
      raise ValueError('Incorrect value for icharg')
    if value > 9: 
      value -= 10
      instance.scf = True
    elif value != -1: instance.scf = False
    self.value = value

  def output_map(self, **kwargs):
    from ..misc import latest_file, copyfile
    from . import files

    icharge = self.value
    if icharge is None: return None
    # some files will be copied.
    if icharge not in [2, 12]: 
      # determines directories to look into.
      vasp = kwargs['vasp']
      outdir = kwargs['outdir']
      hasrestart = getattr(vasp.restart, 'success', False)
      directories = [outdir]
      if hasrestart: directories += [vasp.restart.directory]
      # determines which files exist
      last_wfn = None if icharge in [1, 11]                                    \
                 else latest_file(files.WAVECAR, *directories)
      last_chg = latest_file(files.CHGCAR, *directories) 
      last_pot = None if icharge != 4 else latest_file(files.POT, *directories) 

      # determines icharge depending on file. 
      if last_wfn is not None: icharge = 10 if vasp.nonscf else 0
      elif last_chg is not None: icharge = 11 if vasp.nonscf else 1
      if last_pot is not None and not vasp.nonscf: icharge = 4
      if icharge < 0: return None

      # copies relevant files.
      if last_wfn is not None: copyfile(last_wfn, outdir, nothrow='same')
      if last_chg is not None: copyfile(last_chg, outdir, nothrow='same')
      if last_pot is not None: copyfile(last_pot, outdir, nothrow='same')
    return {self.keyword: icharge}

class IStart(ValueKeyword):
  """ Starting wavefunctions.

      It is best to keep this attribute set to -1, in which case, LaDa takes
      care of copying the relevant files.

        - -1: Automatically determined by LaDA. Depends on the value of restart_
              and the existence of the relevant files.
  
        - 0: Start from scratch.

        - 1: Restart with constant cutoff.
  
        - 2: Restart with constant basis.
  
        - 3: Full restart, including TMPCAR.

      .. note::
      
         Files are copied right before the calculation takes place, not before.

      .. seealso:: ISTART_

      .. _ISTART: http://cms.mpi.univie.ac.at/wiki/index.php/ISTART
      .. _restart: :py:attr:`~lada.vasp.functional.Functional.restart`
  """ 
  keyword = 'ISTART'
  """ VASP keyword """
  aliases = { -1: ['auto', -1], 0: ['scracth', 0],
               1: ['cutoff', 1], 2: ['basis', 2], 3: ['tmpcar', 'full', 3] }
  """ Mapping of aliases. """
  def __init__(self, value=-1): 
    super(IStart, self).__init__(value)

  def output_map(self, **kwargs):
    from ..misc import latest_file, copyfile
    from ..error import RuntimeError
    from . import files

    istart = self.value
    if istart is None: return None
    # some files will be copied.
    if istart != 0:
      # determines directories to look into.
      vasp = kwargs['vasp']
      outdir = kwargs['outdir']
      hasrestart = getattr(vasp.restart, 'success', False)
      directories = [outdir]
      if hasrestart: directories += [vasp.restart.directory]
      # determines which files exist
      last_wfn = latest_file(files.WAVECAR, *directories)
      last_tmp = latest_file(files.TMPCAR, *directories) 

      # determines icharge depending on file. 
      if last_wfn is not None:
        if istart < 0: istart = 1
        else:
          raise RuntimeError( 'Wavefunction file does not exist and ISTART={0}'\
                              .format(istart) )
      if istart == 4 and last_tmp is None:
        raise RuntimeError( 'TMPCAR file does not exist and ISTART={0}'\
                            .format(istart) )
      if last_wfn is not None: copyfile(last_wfn, outdir, nothrow='same')
      if last_tmp is not None: copyfile(last_tmp, outdir, nothrow='same')
    return {self.keyword: istart}

class IStruc(AliasKeyword):
  """ Initial starting structure. """
  aliases = { 'auto': ['auto', 0], 'scratch': ['scracth', 1],
              'contcar': ['contcar', 2], 'restart': ['restart', 3] }
  """ Aliases for the same option. """
  def __init__(self, value='auto'):
    super(IStruc, self).__init__(value=value)
  def output_map(self, **kwargs):
    from ..misc import latest_file, copyfile
    from ..error import RuntimeError, ValueError
    from ..crystal import write
    from . import files

    istruc = self.istruc
    if istruc is None: return None
    # some files will be copied.
    # determines directories to look into.
    vasp = kwargs['vasp']
    outdir = kwargs['outdir']
    hasrestart = getattr(vasp.restart, 'success', False)
    directories = [outdir]
    if hasrestart: directories += [vasp.restart.directory]
    # determines which files exist
    last_poscar = latest_file(files.POSCAR, *directories)
    last_contcar = latest_file(files.CONTCAR, *directories) 
    if last_contcar is not None and last_contcar is not None:
      earliest = stat(last_contcar).st_mtime > stat(last_poscar).st_mtime
      earliest = last_contcar if earliest else last_poscar
    else: earliest = last_contcar if last_contcar is not None else last_poscar

    if istruc == 'scratch' or (last_poscar is None and last_contcar is None):
      structure = kwargs['structure']
      if len(structure) == 0: raise ValueError('Structure is empty')
      if structure.scale < 1e-8: raise ValueError('Structure scale is zero')
      if structure.volume < 1e-8: raise ValueError('Structure volume is zero')
      write.poscar(structure)
    if istruc == 'contcar': 
      if not last_contcar: raise ValueError('No contcar file to copy')
      copyfile(last_contcar, '.', nothrow='same')
    elif istruc == 'restart':
      if not last_
      
      raise ValueError(
    # determines icharge depending on file. 
      if last_wfn is not None:
        if istart < 0: istart = 1
        else:
          raise RuntimeError( 'Wavefunction file does not exist and ISTART={0}'\
                              .format(istart) )
      if istart == 4 and last_tmp is None:
        raise RuntimeError( 'TMPCAR file does not exist and ISTART={0}'\
                            .format(istart) )
      if last_wfn is not None: copyfile(last_wfn, outdir, nothrow='same')
      if last_tmp is not None: copyfile(last_tmp, outdir, nothrow='same')
    return {self.keyword: istart}



class PartialRestart(ValueKeyword):
  """ Restart from previous run.
      
      It is either an vasp extraction object of some kind, or None.  In the
      latter case, the calculation starts from scratch.  However, if an
      extraction object exists *and* the calculation it refers to was
      successfull, then it will check whether WAVECAR and CHGCAR exist. It also
      checks whether :py:attr:`nonscf <incar.Incar.nonscf>` is True or False,
      and sets ICHARG_ accordingly. The CONTCAR file is *never* copied from the
      previous run. For an alternate behavior, see :py:class:`Restart`.

      .. seealso:: ICHARG_, ISTART_, :py:class:`Restart`
  """
  def __init__(self, value): super(PartialRestart, self).__init__(value=value)

  def output_map(self, **kwargs):
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
      copyfile(join(self.value.directory, files.WAVEDER), files.WAVEDER,
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),
               nocopyempty=True) 
      copyfile(join(self.value.directory, files.TMPCAR), files.TMPCAR,
               nothrow='same exists', symlink=getattr(kwargs["vasp"], 'symlink', False),
               nocopyempty=True) 
      if kwargs['vasp'].lsorbit == True: kwargs['vasp'].nbands = 2*self.value.nbands 
    return None
