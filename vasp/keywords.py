from ..functools.keyword import BoolKeyword as BaseBoolKeyword, ValueKeyword,  \
                                TypedKeyword
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
