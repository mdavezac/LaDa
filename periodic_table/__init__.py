""" Physical quantities of elements. 

    Atomic quantities can be used as in:

    .. python::
    
      import lada.periodic_table
      periodic_table.Au.atomic_weight
      periodic_table.Gold.atomic_weight

    Available quantities can be found in `lada.periodic_table.Atom`. Some of
    these are set to None when either meaningless or not available.
"""
__docformat__ = "restructuredtext en"
import quantities as qt


class Element(object):
  """ Contains atomic data for single element. 
  
      Data is taken from the `webelements`_ website.

      .. _webelements: http://www.webelements.com
  """
  def __init__(self, **kwargs):
    """ Initializes atoms """
    self.symbol                 = kwargs.pop('symbol', None)
    """ Atomic symbol. """
    self.name                   = kwargs.pop('name', None)
    """ Name of the element. """
    self.atomic_weight          = kwargs.pop('atomic_weight', None)
    """ Atomic weight (dimensionless) of the element. """
    self.atomic_number          = kwargs.pop('atomic_number', None)
    """ Atomic number (dimensionless) of the element. """
    self.pauling                = kwargs.pop('pauling', None)
    """ Pauling electronegativity (Pauling scale) of the element. """
    self.sanderson              = kwargs.pop('sanderson', None)
    """ Sanderson electronegativity (Pauling scale) of the element. """
    self.allred_rochow          = kwargs.pop('allred_rochow', None)
    """ Allred-Rochow electronegativity (Pauling scale) of the element. """
    self.mulliken_jaffe         = kwargs.pop('mulliken_jaffe', None)
    """ Mulliken-Jaffe electronegativity (Pauling scale) of the element. """
    self.allen                  = kwargs.pop('allen', None)
    """ Allen electronegativity (Pauling scale) of the element. """
    self.electron_affinity      = kwargs.pop('electron_affinity', None)
    """ Electron affinity of the element (kJ per mol). """
    self.ionization_energies    = kwargs.pop('ionization_energies', None)
    """ Known Ionization energies of the element (kJ per mol).
    
        All ionization energies known to www.webelements.com are listed, from
        first to last.
    """
    self.atomic_radius          = kwargs.pop('atomic_radius', None)
    """ Empirical atomic radius. """
    self.covalent_radius        = kwargs.pop('covalent_radius', None)
    """ Covalent bond radius. """
    self.single_bond_radius     = kwargs.pop('single_bond_radius', None)
    """ Single covalent-bond  radius. """
    self.double_bond_radius     = kwargs.pop('double_bond_radius', None)
    """ Double covalent-bond  radius. """
    self.triple_bond_radius     = kwargs.pop('triple_bond_radius', None)
    """ Triple covalent-bond  radius. """
    self.van_der_waals_radius   = kwargs.pop('van_der_waals_radius', None)
    """ van der Walls radius. """
    self.fusion                 = kwargs.pop('fusion', None)
    """ Enthalpy of fusion. """
    self.vaporization           = kwargs.pop('vaporization', None)
    """ Enthalpy of vaporization. """
    self.atomization            = kwargs.pop('atomization', None)
    """ Enthalpy of atomization. """
    self.melting_point          = kwargs.pop('melting_point', None)
    """ Melting point of the elemental solid. """
    self.boiling_point          = kwargs.pop('boiling_point', None)
    """ Boiling point of the elemental solid. """
    self.critical_temperature   = kwargs.pop('critical_temperature', None)
    """ Critical temperature of the elemental solid. """
    self.thermal_conductivity   = kwargs.pop('thermal_conductivity', None)
    """ Thermal conductivity of the elemental solid. """
    self.thermal_expansion      = kwargs.pop('thermal_expansion', None)
    """ Coefficient of lineary thermal expansion of the elemental solid at 273K. """
    self.density                = kwargs.pop('density', None)
    """ Density of the elemental solid. """
    self.molar_volume           = kwargs.pop('molar_volume', None)
    """ Molar volume of the element at 298K. """
    self.sound_velocity         = kwargs.pop('sound_velocity', None)
    """ Velocity of sound in the element at 298K. """
    self.young_modulus          = kwargs.pop('young_modulus', None)
    """ Young modulus of the elemental solid. """
    self.rigidity_modulus       = kwargs.pop('rigidity_modulus', None)
    """ Rigidity modulus of the elemental solid. """
    self.bulk_modulus           = kwargs.pop('bulk_modulus', None)
    """ Bulk modulus ratio of the elemental solid. """
    self.poisson_ratio          = kwargs.pop('poisson_ratio', None)
    """ Poisson ratio of the elemental solid. """
    self.electrical_resistivity = kwargs.pop('electical_resistivity', None)
    """ Electrical Resistivity ratio of the elemental solid. """
    self.pettifor               = kwargs.pop('pettifor', None)
    """ This element on the `Pettifor` scale.

        `Pettifor`_'s is an artificial scale designed to parameterizes
        a two-dimensional structure map of binary AB compounds. 
  
        References
        ==========
          .. _Pettifor : D.G. Pettifor, Solid. Stat. Comm., *51* 31-34 (1984).
    """
    self.orbital_radii        = kwargs.pop('orbital_radii', None)
    """ `Orbital`_ radii of this element.

        The orbital radii can be defined for s, p, and d orbitals using
        psudo-potential wavefunctions. 
        
        References
        ==========
       
        .. _Orbital radii : Alex Zunger, PRB *22* 5839-5872 (1980),
            http://dx.doi.org/10.1103/PhysRevB.22.5839
    """



  def __str__(self):
    return "{0} ({1}): n={2}, Z={3}"\
           .format(self.name, self.symbol, self.atomic_number, self.atomic_weight)
  def __repr__(self):
    result = {}
    for name in self.__dict__:
      if name[0] == '_': continue
      s = getattr(self, name)
      if s != None: result[name] = s
    return "{0}(**{1})".format(self.__class__.__name__, result)



def iterate():
  """ Iterates through all elements. """
  for name in _elements.symbols:
    yield globals()[name] 


from _create_data import *
import _elements
__all__ = list(_elements.symbols)
__all__.extend(['Element', 'iterate'])
locals().update(_elements.elements)
