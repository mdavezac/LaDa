""" Physical quantities of elements. 

    Atomic quantities can be used as in:
      >>> import lada.periodic_table
      >>> periodic_table.Au.atomic_weight
      >>> periodic_table.Gold.atomic_weight
    Available quantities can be found in `lada.periodic_table.Atom`. Some of
    these are set to None when either meaningless or not available.
"""
__docformat__ = "restructuredtext en"
import quantities as qt
# from _elements import *

# __dir__ = _elements.__dir__ 
# __dir__.extend(['Element'])

class Element(object):
  """ Contains atomic data for single element. 
  
      Data is taken from www.webelements.com.
  """
  __slots__ = ['symbol', 'name', 'atomic_weight', 'atomic_number', 'pauling', 'sanderson', \
               'allred_rochow', 'mulliken_jaffe', 'allen', 'electron_affinity', \
               'ionization_energies', 'atomic_radius', 'covalent_radius', 'single_bond_radius', \
               'double_bond_radius', 'triple_bond_radius', 'van_der_waals_radius', \
               'fusion', 'vaporization', 'atomization', 'melting_point', 'boiling_point', \
               'critical_temperature', 'thermal_conductivity', 'thermal_expansion', 'density',\
               'molar_volume', 'sound_velocity', 'young_modulus', 'rigidity_modulus',\
               'bulk_modulus', 'poisson_ratio', 'electrical_resistivity' ]
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



  def __str__(self):
    return "{0} ({1}): n={2}, Z={3}"\
           .format(self.name, self.symbol, self.atomic_number, self.atomic_weight)
  def __repr__(self):
    result = {}
    for name in self.__slots__:
      s = getattr(self, name)
      if s != None: result[name] = s
    return "{0} = {1}({2})".format(self.name, self.__class__.__name__, result)


def _download_files():
  """ Downloads data from webelements.com. """
  import urllib
  from os import makedirs
  from os.path import exists, join
  
  atom_list = ['Ruthenium', 'Rhenium', 'Rutherfordium', 'Radium', 'Rubidium',
    'Radon', 'Rhodium', 'Beryllium', 'Barium', 'Bohrium', 'Bismuth',
    'Berkelium', 'Bromine', 'Hydrogen', 'Phosphorus', 'Osmium', 'Mercury',
    'Germanium', 'Gadolinium', 'Gallium', 'Ununbium', 'Praseodymium',
    'Platinum', 'Plutonium', 'Carbon', 'Lead', 'Protactinium', 'Palladium',
    'Xenon', 'Polonium', 'Promethium', 'Hassium',
    'Holmium', 'Hafnium', 'Molybdenum', 'Helium', 'Mendelevium', 'Magnesium',
    'Potassium', 'Manganese', 'Oxygen', 'Meitnerium', 'Sulfur', 'Tungsten',
    'Zinc', 'Europium', 'Einsteinium', 'Erbium', 'Nickel', 'Nobelium',
    'Sodium', 'Niobium', 'Neodymium', 'Neon', 'Neptunium', 'Francium', 'Iron',
    'Fermium', 'Boron', 'Fluorine', 'Strontium', 'Nitrogen', 'Krypton',
    'Silicon', 'Tin', 'Samarium', 'Vanadium', 'Scandium', 'Antimony',
    'Seaborgium', 'Selenium', 'Cobalt', 'Curium', 'Chlorine', 'Calcium',
    'Californium', 'Cerium', 'Cadmium', 'Thulium', 'Caesium', 'Chromium',
    'Copper', 'Lanthanum', 'Lithium', 'Thallium', 'Lutetium', 'Lawrencium',
    'Thorium', 'Titanium', 'Tellurium', 'Terbium', 'Technetium', 'Tantalum',
    'Ytterbium', 'Dubnium', 'Zirconium', 'Dysprosium', 'Iodine', 'Uranium',
    'Yttrium', 'Actinium', 'Silver', 'Iridium', 'Americium', 'Aluminium',
    'Arsenic', 'Argon', 'Gold', 'Astatine', 'Indium', 'Darmstadtium', 'Copernicium']

  if not exists("elements"): makedirs("elements")
  for name in atom_list: 
    file = urllib.urlopen("http://www.webelements.com/{0}".format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name), "w") as out: out.write(string)
    file = urllib.urlopen("http://www.webelements.com/{0}/atoms.html".format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name + "_atoms.html"), "w") as out: out.write(string)
    file = urllib.urlopen( "http://www.webelements.com/{0}/electronegativity.html"\
                           .format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name + "_electronegativity.html"), "w") as out: out.write(string)
    file = urllib.urlopen( "http://www.webelements.com/{0}/atom_sizes.html"\
                           .format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name + "_atom_sizes.html"), "w") as out: out.write(string)
    file = urllib.urlopen( "http://www.webelements.com/{0}/thermochemistry.html"\
                           .format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name + "_thermochemistry.html"), "w") as out: out.write(string)
    file = urllib.urlopen( "http://www.webelements.com/{0}/physics.html"\
                           .format(name.lower()))
    string = file.read()
    file.close()
    with open(join("elements", name + "_physics.html"), "w") as out: out.write(string)


def _create_elements_py(filename="_elements.py"):
  """ Gets data from webelements.com and creates _elements.py. """
  import re
  import urllib
  from os.path import exists, join
  from BeautifulSoup import BeautifulSoup, HTMLParseError
  import quantities as pq

  atom_list = ['Ruthenium', 'Rhenium', 'Rutherfordium', 'Radium', 'Rubidium',
    'Radon', 'Rhodium', 'Beryllium', 'Barium', 'Bohrium', 'Bismuth',
    'Berkelium', 'Bromine', 'Hydrogen', 'Phosphorus', 'Osmium', 'Mercury',
    'Germanium', 'Gadolinium', 'Gallium', 'Ununbium', 'Praseodymium',
    'Platinum', 'Plutonium', 'Carbon', 'Lead', 'Protactinium', 'Palladium',
    'Xenon', 'Polonium', 'Promethium', 'Hassium',
    'Holmium', 'Hafnium', 'Molybdenum', 'Helium', 'Mendelevium', 'Magnesium',
    'Potassium', 'Manganese', 'Oxygen', 'Meitnerium', 'Sulfur', 'Tungsten',
    'Zinc', 'Europium', 'Einsteinium', 'Erbium', 'Nickel', 'Nobelium',
    'Sodium', 'Niobium', 'Neodymium', 'Neon', 'Neptunium', 'Francium', 'Iron',
    'Fermium', 'Boron', 'Fluorine', 'Strontium', 'Nitrogen', 'Krypton',
    'Silicon', 'Tin', 'Samarium', 'Vanadium', 'Scandium', 'Antimony',
    'Seaborgium', 'Selenium', 'Cobalt', 'Curium', 'Chlorine', 'Calcium',
    'Californium', 'Cerium', 'Cadmium', 'Thulium', 'Caesium', 'Chromium',
    'Copper', 'Lanthanum', 'Lithium', 'Thallium', 'Lutetium', 'Lawrencium',
    'Thorium', 'Titanium', 'Tellurium', 'Terbium', 'Technetium', 'Tantalum',
    'Ytterbium', 'Dubnium', 'Zirconium', 'Dysprosium', 'Iodine', 'Uranium',
    'Yttrium', 'Actinium', 'Silver', 'Iridium', 'Americium', 'Aluminium',
    'Arsenic', 'Argon', 'Gold', 'Astatine', 'Indium']

  re_swf = re.compile("(rainbow|NI3|volcano|\_flash|K\_H2O).swf\s*(?!\")")
  re_atomweight = re.compile(":\s*\[?\s*(\d+(?:\.\d+)?)\s*\]?")
  results = {}
  for name in atom_list: 

    # first opens and reads file.
    if not exists(join("elements", name)): 
      file = urllib.urlopen("http://www.webelements.com/{0}".format(name.lower()))
      string = file.read()
      file.close()
    else:
      with open(join("elements", name), "r") as file: string = file.read()
    string = string.replace("alt\"", "alt=\"")
    soup = BeautifulSoup(re.sub(re_swf,"rainbow.swf\"",string))

    atom = Element(name=name)
    atom.symbol = soup.findChild( name="a", attrs={"title": "Element names and symbols"},\
                                  text=" Symbol").parent.parent.contents[1].split()[1]
    atom.atomic_number = soup.findChild(name="a", attrs={"title": "Element atomic numbers"})\
                                       .parent.contents[-1].split()[1]
    atom.atomic_number = int(atom.atomic_number)
    atom.atomic_weight = soup.findChild(name="a", attrs={"title": "Element atomic weights"})\
                                       .parent.prettify()
    found = re_atomweight.search(atom.atomic_weight)
    if found == None: print name
    else: atom.atomic_weight = float(found.group(1))

    
    # ionization stuff
    if not exists(join("elements", name + "_atoms.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/atoms.html".format(name.lower()))
      string = file.read()
      file.close()
    else: 
      with open(join("elements", name + "_atoms.html"), "r") as file: string = file.read()
    soup = BeautifulSoup(string) 
    # electron affinity
    found = re.search("of\s+{0}\s+is\s+(\S+)".format(name.lower()), string)
    if found.group(1) == "no": atom.electron_affinity = None
    else: atom.electron_affinity = float(found.group(1)) * pq.kilo * pq.J / pq.mol
    # ionization energies
    energies = []
    for child in soup.findChild(name="table", attrs={"class":"chemistry-data"})\
                     .findChildren(name='td'):
      energies.append(float(child.string) * pq.kilo * pq.J / pq.mol)
    atom.ionization_energies = energies if len(energies) > 0 else None


    # electronegativities.
    if not exists(join("elements", name + "_electronegativity.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/electronegativity.html"\
                            .format(name.lower()))
      string = file.read()
      file.close()
    else: 
      with open(join("elements", name + "_electronegativity.html"), "r") as file:
          string = file.read()
    soup = BeautifulSoup(string) 
    attrs = { "href": "../periodicity/electronegativity_pauling/",\
              "title": "View definition and pictures showing periodicity "\
                       "of Pauling electronegativity"}
    pauling = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
    pauling = pauling.split()[0]
    atom.pauling = float(pauling) if pauling != "no" else None

    attrs = { "href": "../periodicity/electronegativity_sanderson/" }
    sanderson = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
    sanderson = sanderson.split()[0]
    atom.sanderson = float(sanderson) if sanderson != "no" else None

    attrs = { "href": "../periodicity/electroneg_allred_rochow/" }
    allred_rochow = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
    allred_rochow = allred_rochow.split()[0]
    atom.allred_rochow = float(allred_rochow) if allred_rochow != "no" else None

    attrs = { "href": "../periodicity/electroneg_mulliken_jaffe/" }
    mulliken_jaffe = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1]
    if name in ["Germanium", "Gallium", "Carbon", "Lead", "Boron", "Silicon", "Tin",\
                "Thallium", "Aluminium", "Indium"]: 
      mulliken_jaffe = mulliken_jaffe.contents[0]
    else: mulliken_jaffe = mulliken_jaffe.string
    mulliken_jaffe = mulliken_jaffe.split()[0]
    atom.mulliken_jaffe = float(mulliken_jaffe) if mulliken_jaffe != "no" else None

    attrs = { "href": "../periodicity/electronegativity_allen/" }
    allen = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
    allen = allen.split()[0]
    atom.allen = float(allen) if allen != "no" else None
    
    # atom sizes
    if not exists(join("elements", name + "_atom_sizes.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/atom_sizes.html"\
                            .format(name.lower()))
      string = file.read()
      file.close()
    else: 
      with open(join("elements", name + "_atom_sizes.html"), "r") as file:
          string = file.read()
    soup = BeautifulSoup(string) 
    
    # atomic radius
    attrs = { "href": "../periodicity/atomic_radius_empirical/" }
    atomic_radius = soup.findChild(name="a", attrs=attrs).parent.contents[-1].split()[1]
    if atomic_radius != "no":
      atom.atomic_radius = float(atomic_radius) * pq.picometre 
    
    attrs = { "href": "../periodicity/covalent_radius_2008/" }
    covalent_radius = soup.findChild(name="a", attrs=attrs).parent.contents[-1].split()[1]
    atom.covalent_radius = float(covalent_radius) * pq.picometre if covalent_radius != "no" else None

    attrs = { "href": "../periodicity/radii_covalent_single/" }
    single_bond_radius = soup.findChild(name="a", attrs=attrs)
    if single_bond_radius != None:
      single_bond_radius = single_bond_radius.parent.contents[-1].split()[1]
      if single_bond_radius != "no": 
        atom.single_bond_radius = float(single_bond_radius) * pq.picometre

    attrs = { "href": "../periodicity/radii_covalent_double/" }
    double_bond_radius = soup.findChild(name="a", attrs=attrs)
    if double_bond_radius != None:
      double_bond_radius = double_bond_radius.parent.contents[-1].split()[1]
      if double_bond_radius != "no": 
        atom.double_bond_radius = float(double_bond_radius) * pq.picometre

    attrs = { "href": "../periodicity/radii_covalent_triple/" }
    triple_bond_radius = soup.findChild(name="a", attrs=attrs)
    if triple_bond_radius != None:
      triple_bond_radius = triple_bond_radius.parent.contents[-1].split()[1]
      if triple_bond_radius != "no": 
        atom.triple_bond_radius = float(triple_bond_radius) * pq.picometre

    attrs = { "href": "../periodicity/van_der_waals_radius/" }
    van_der_waals_radius = soup.findChild(name="a", attrs=attrs)
    if van_der_waals_radius != None:
      van_der_waals_radius = van_der_waals_radius.parent.contents[-1].split()[1]
      if van_der_waals_radius != "no": 
        atom.van_der_waals_radius = float(van_der_waals_radius) * pq.picometre

    # thermochemistry
    if not exists(join("elements", name + "_thermochemistry.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/thermochemistry.html"\
                            .format(name.lower()))
      string = file.read()
      file.close()
    else: 
      with open(join("elements", name + "_thermochemistry.html"), "r") as file:
          string = file.read()
    soup = BeautifulSoup(string) 
    
    attrs = { "href": "../periodicity/enthalpy_fusion/" }
    fusion = soup.findChild(name="a", attrs=attrs).parent.prettify()
    fusion = re.search(":\s*(?:about)?\s*(\S+)", fusion)
    if fusion != None and fusion.group(1) != "no":
      atom.fusion = float(fusion.group(1)) * pq.kilo * pq.J / pq.mol 

    attrs = { "href": "../periodicity/enthalpy_vaporisation/" }
    vaporization = soup.findChild(name="a", attrs=attrs).parent.prettify()
    vaporization = re.search(":\s*(?:about)?\s*(\S+)", vaporization)
    if vaporization != None and vaporization.group(1) != "no":
      atom.vaporization = float(vaporization.group(1)) * pq.kilo * pq.J / pq.mol 

    attrs = { "href": "../periodicity/enthalpy_atomisation/" }
    atomization = soup.findChild(name="a", attrs=attrs).parent.prettify()
    atomization = re.search(":\s*(?:about)?\s*(\S+)", atomization)
    if atomization != None and atomization.group(1) != "no":
      atom.atomization = float(atomization.group(1)) * pq.kilo * pq.J / pq.mol 

    # physics
    if not exists(join("elements", name + "_physics.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/physics.html"\
                            .format(name.lower()))
      string = file.read()
      file.close()
    else: 
      with open(join("elements", name + "_physics.html"), "r") as file:
          string = file.read()
    soup = BeautifulSoup(string) 

    attrs = { "href": "../periodicity/melting_point/" }
    melting_point = soup.findChild(name="a", attrs=attrs).parent.prettify()
    melting_point = re.search(":\s*(?:\(white P\)|about|maybe about)?\s*(\S+)", melting_point)
    if melting_point != None and melting_point.group(1) != "no":
      atom.melting_point = float(melting_point.group(1)) * pq.Kelvin

    attrs = { "href": "../periodicity/boiling_point/" }
    boiling_point = soup.findChild(name="a", attrs=attrs).parent.prettify()
    boiling_point = re.search(":\s*(?:about)?\s*(\S+)", boiling_point)
    if boiling_point != None and boiling_point.group(1) != "no":
      atom.boiling_point = float(boiling_point.group(1)) * pq.Kelvin

    attrs = { "href": "../periodicity/critical_temperature/" }
    critical_temperature = soup.findChild(name="a", attrs=attrs).parent.prettify()
    critical_temperature = re.search(":\s*(?:about)?\s*(\S+)", critical_temperature)
    if critical_temperature != None and critical_temperature.group(1) != "no":
      atom.critical_temperature = float(critical_temperature.group(1)) * pq.Kelvin

    attrs = { "href": "../periodicity/thermal_conductivity/" }
    thermal_conductivity = soup.findChild(name="a", attrs=attrs).parent.prettify()
    thermal_conductivity = re.search(":\s*(?:about)?\s*(\S+)", thermal_conductivity)
    if thermal_conductivity != None and thermal_conductivity.group(1) != "no":
      atom.thermal_conductivity = float(thermal_conductivity.group(1)) * pq.W / pq.m / pq.K

    attrs = { "href": "../periodicity/coeff_thermal_expansion/" }
    thermal_expansion = soup.findChild(name="a", attrs=attrs).parent.prettify()
    thermal_expansion = re.search(":\s*(?:about)?\s*(\S+)", thermal_expansion)
    if thermal_expansion != None and thermal_expansion.group(1) != "no":
      atom.thermal_expansion = float(thermal_expansion.group(1)) * pq.micro / pq.K

    attrs = { "href": "../periodicity/density/" }
    density = soup.findChild(name="a", attrs=attrs).parent.prettify()
    density = re.search(":\s*(?:about)?\s*(\S+)", density)
    if density != None and density.group(1) != "no":
      atom.density = float(density.group(1)) / 1000 * pq.g * pq.cm**3

    attrs = { "href": "../periodicity/molar_volume/" }
    molar_volume = soup.findChild(name="a", attrs=attrs).parent.prettify()
    molar_volume = re.search(":\s*(?:about)?\s*(\S+)", molar_volume)
    if molar_volume != None and molar_volume.group(1) != "no":
      atom.molar_volume = float(molar_volume.group(1)) * pq.cm**3 / pq.mol

    attrs = { "href": "../periodicity/velocity_sound/" }
    sound_velocity = soup.findChild(name="a", attrs=attrs).parent.prettify()
    sound_velocity = re.search(":\s*(?:about)?\s*(\S+)", sound_velocity)
    if sound_velocity != None and sound_velocity.group(1) != "no":
      atom.sound_velocity = float(sound_velocity.group(1)) * pq.m / pq.s

    attrs = { "href": "../periodicity/youngs_modulus/" }
    young_modulus = soup.findChild(name="a", attrs=attrs).parent.prettify()
    young_modulus = re.search(":\s*(?:about)?\s*(\S+)", young_modulus)
    if young_modulus != None and young_modulus.group(1) != "no":
      atom.young_modulus = float(young_modulus.group(1)) * pq.GPa

    attrs = { "href": "../periodicity/rigidity_modulus/" }
    rigidity_modulus = soup.findChild(name="a", attrs=attrs).parent.prettify()
    rigidity_modulus = re.search(":\s*(?:about)?\s*(\S+)", rigidity_modulus)
    if rigidity_modulus != None and rigidity_modulus.group(1) != "no":
      atom.rigidity_modulus = float(rigidity_modulus.group(1)) * pq.GPa
    
    attrs = { "href": "../periodicity/bulk_modulus/" }
    bulk_modulus = soup.findChild(name="a", attrs=attrs).parent.prettify()
    bulk_modulus = re.search(":\s*(?:about)?\s*(\S+)", bulk_modulus)
    if bulk_modulus != None and bulk_modulus.group(1) != "no":
      atom.bulk_modulus = float(bulk_modulus.group(1)) * pq.GPa
    
    attrs = { "href": "../periodicity/poissons_ratio/" }
    poisson_ratio = soup.findChild(name="a", attrs=attrs).parent.prettify()
    poisson_ratio = re.search(":\s*(?:about)?\s*(\S+)", poisson_ratio)
    if poisson_ratio != None and poisson_ratio.group(1) != "no":
      atom.poisson_ratio = float(poisson_ratio.group(1)) * pq.dimensionless
    
    attrs = { "href": "../periodicity/electrical_resistivity/" }
    electrical_resistivity = soup.findChild(name="a", attrs=attrs).parent.prettify()
    electrical_resistivity = re.search(":\s*(?:about)?\s*(\d+(?:\.\d+)?)", electrical_resistivity)
    if electrical_resistivity != None and electrical_resistivity.group(1) not in ["no", "&gt;"]:
      atom.electrical_resistivity = float(electrical_resistivity.group(1)) * 1e-8 * pq.ohm * pq.m
    print repr(atom)
