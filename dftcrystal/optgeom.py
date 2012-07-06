from input import AttrBlock, BoolKeyword, TypedKeyword

class GeometryOpt(BoolKeyword):
  """ Geometry optimization keyword. 

      When set to True, makes sure that:

        1. MAXCYCLE is non-zero. If it is, sets it to 100.
        2. only one of FULLOPTG, ITATOCELL, ITREDUN, CELLONLY is True.
        3. CVOLOPT is False unless this is a FULLOPTG or CELLONLY calculation.

      (2) and (3) are decided by the keyword argument in input.
  """
  def __init__(self, keyword, value=False):
    """ Initializes a geometry optimization keyword. """
    self.keyword = keyword
    """ CRYSTAL keyword. """
    super(GeometryOpt, self).__init__(value=value)

  def __set__(self, instance, value):
    """ Sets keyword to appear or not. 

        Also clears attendant conditions described in :py:class:`GeometryOpt`.
    """
    super(GeometryOpt, self).__set__(instance, value)
    if self.value:
      if instance.maxcycle is not None and instance.maxcycle < 1:
        instance.maxcycle = 100
      if self.keyword != 'fulloptg':  instance.fulloptg  = False
      if self.keyword != 'cellonly':  instance.cellonly  = False
      if self.keyword != 'itatocell': instance.itatocell = False
      if self.keyword != 'interdun':  instance.interdun  = False
      if self.keyword != 'fulloptg' and self.keyword != 'cellonly':
        instance.cvolopt = False

class CVolOpt(BoolKeyword):
  """ Constant volume optimization keyword. 
  
      Only appears if FULLOPTG or CELLONLY exist.
  """
  keyword = 'cvolopt'
  """ Crystal input keyword """
  def print_input(self, **kwargs):
    """ Only prints if FULLOPTG or CELLONLY exist. """
    optgeom = kwargs['crystal'].optgeom
    if optgeom.fulloptg == False and optgeom.cellonly == False:
      return None
    return super(CVolOpt, self).print_input(**kwargs)

class MaxCycle(TypedKeyword):
  """ Maxcycle input.

      It is expected to be an integer or None.
      If changed from None or 0 to an integer greater than zero, and no
      optimization method has yet been selected, set the optimization method to
      fulloptg.
  """
  keyword = 'maxcycle'
  """ CRYSTAL keyword. """
  def __init__(self, value=None):
    """ Creates MAXCYCLE keyword. """
    super(MaxCycle, self).__init__(keyword=None, value=value, type=int)
  def __set__(self, instance, value):
    super(MaxCycle, self).__set__(instance, value)
    if self._value is not None and self._value > 0 and instance.static:
      instance.fulloptg = True

class ExclAttrBlock(AttrBlock):
  """ An attribute block set up to exclude others. 
  
      Expects both "keyword" and "excludegroup" attributes (or class
      attributes) to exist in derived instances.
  """
  excludegroup = 'optgeom', 'freqcalc', 'anharm', 'confcnt', 'cphf',           \
                 'elastcon', 'eos'
  """ Groups from which only one should be enabled. """
  def __init__(self):
    """ Initializes the exclusive attribute block. """
    super(ExclAttrBlock, self).__init__()
    self.enabled = False
    self._parent = None
    """ Weak-reference to parent instance. """

  @property
  def enabled(self):
    """ True  if this block is enabled. """
    return self._enabled
  @enabled.setter
  def enabled(self, value):
    self._enabled = value == True
    # disable other instances. 
    if self._enabled and self._parent is not None                               \
       and len(getattr(self, 'excludegroup', ())) > 0:
      parent = self._parent()
      for u in self.excludegroup:
        if u == self.keyword: continue
        if not hasattr(parent, u): continue
        inst = getattr(parent, u)
        if not hasattr(inst, 'enabled'): continue
        inst.enabled = False
  def __get__(self, instance, owner=None):
    """ Sets up weak ref to parent. """
    from weakref import ref
    self._parent = ref(instance)
    return self

  def read_input(self, tree, owner=None):
    """ Reads from input. """
    if owner is not None:
      from weakref import ref
      self._parent = ref(owner)
    self.enabled = True
    return super(ExclAttrBlock, self).read_input(tree, owner)

  def print_input(self, **kwargs):
    """ Does not print if static. """
    if not self.enabled: return None
    return super(ExclAttrBlock, self).print_input(**kwargs)

  def __getstate__(self):
    d = self.__dict__.copy()
    d['_parent'] = None
    return d
  def __setstate__(self, value):
    self.__dict__.update(value)

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly output. """
    # make sure enabled is printed.
    if defaults is not None: defaults.enabled = not self.enabled
    return super(ExclAttrBlock, self).__ui_repr__( imports, name,              \
                                                   defaults, exclude )

class OptGeom(ExclAttrBlock):
  """ Geometry Optimization block. """
  keyword = 'optgeom'
  """ CRYSTAL input keyword """
  def __init__(self): 
    """ Creates an optimization block. """
    from .input import QuantityKeyword
    from quantities import UnitQuantity, hartree, angstrom
    super(OptGeom, self).__init__()
 
    self.maxcycle   = MaxCycle()
    """ Maxium number of iterations in geometry optimization loop.
    
        It is expected to be an integer or None.
        If changed from None or 0 to an integer greater than zero, and no
        optimization method has yet been selected, set the optimization method
        to fulloptg.
    """
    self.fulloptg   = GeometryOpt('fulloptg')
    """ Full optimization """
    self.cellonly   = GeometryOpt('cellonly')
    """ Cell-shape optimization only """
    self.itatocell  = GeometryOpt('itatocell')
    """ Iterative cell-shape/atom optimization """
    self.interdun   = GeometryOpt('interdun')
    """ Constrained optimization """
    self.cvolopt    = CVolOpt()
    """ Constant volume optimization """
    self.toldee     = TypedKeyword('toldee', int)
    self.toldeg     = TypedKeyword('toldeg', float)
    self.toldex     = TypedKeyword('toldex', float)
    bohr = UnitQuantity('crystal_bhor', 0.5291772083*angstrom, symbol='bhor')
    """ Bhor radius as defined by CRYSTAL """
    self.extpress   = QuantityKeyword(units=hartree/bohr**3)
    """ Hydrostatic pressure in hartree*bohr^-3 """
    self.extstress  = QuantityKeyword(units=hartree/bohr**3, shape=(3,3))
    """ External stress in hartree*bhor^-3 """

  @property
  def static(self): 
    """ Sets to static calculation.
    
        Calculations are static either if maxcycle is 0 or if no optimization
        method is given.

        Calculations are set static - if thery already are not - by setting maxcycle to 0.
        Calculations are made non-static by setting maxcycle to None (e.g.
        CRYSTAL default) if maxcycle is 0, and by setting the optimization
        method to fulloptg if none are yet selected.
    """
    return (self.maxcycle is not None and self.maxcycle < 1)                   \
           or ( self.fulloptg == False and self.cellonly == False              \
                and self.itatocell == False and self.interdun == False )
  @static.setter
  def static(self, value):
    """ If True, makes calculation static. """
    if value == True:
      if self.static: return
      self.maxcycle = 0
    else:
      # Makes sure we are doing some cycles.
      if self.maxcycle is None or self.maxcycle < 1: self.maxcyle = None
      # Makes sure some optimization method is given.
      if self.fulloptg == False and self.cellonly == False                     \
         and self.itatocell == False and self.interdun == False: 
        self.fullopg = True

  def __ui_repr__(self, imports, name=None, defaults=None, exclude=None):
    """ User-friendly output. """
    return super(OptGeom, self).__ui_repr__( imports, name, 
                                             defaults, ['static'] )
  

  def print_input(self, **kwargs):
    """ Does not print if static. """
    if self.static: return None
    return super(OptGeom, self).print_input(**kwargs)
