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
               'boiling_point', 'melting_point', \
               'density']
  def __init__(self, **kwargs):
    """ Initializes atoms """
    self.symbol              = kwargs.pop('symbol', None)
    """ Atomic symbol. """
    self.name                = kwargs.pop('nam', None)
    """ Name of the element. """
    self.atomic_weight       = kwargs.pop('atomic_weight', None)
    """ Atomic weight (dimensionless) of the element. """
    self.atomic_number       = kwargs.pop('atomic_number', None)
    """ Atomic number (dimensionless) of the element. """
    self.pauling             = kwargs.pop('pauling', None)
    """ Pauling electronegativity (Pauling scale) of the element. """
    self.sanderson           = kwargs.pop('sanderson', None)
    """ Sanderson electronegativity (Pauling scale) of the element. """
    self.allred_rochow       = kwargs.pop('allred_rochow', None)
    """ Allred-Rochow electronegativity (Pauling scale) of the element. """
    self.mulliken_jaffe      = kwargs.pop('mulliken_jaffe', None)
    """ Mulliken-Jaffe electronegativity (Pauling scale) of the element. """
    self.allen               = kwargs.pop('allen', None)
    """ Allen electronegativity (Pauling scale) of the element. """
    self.electron_affinity   = kwargs.pop('electron_affinity', None)
    """ Electron affinity of the element (kJ per mol). """
    self.ionization_energies = kwargs.pop('ionization_energies', None)
    """ Known Ionization energies of the element (kJ per mol).
    
        All ionization energies known to www.webelements.com are listed, from
        first to last.
    """
    self.atomic_radius       = kwargs.pop('atomic_radius', None)
    """ Empirical atomic radius. """
    self.covalent_radius       = kwargs.pop('covalent_radius', None)
    """ Covalent bond radius. """
    self.single_bond_radius       = kwargs.pop('single_bond_radius', None)
    """ Single covalent-bond  radius. """
    self.double_bond_radius       = kwargs.pop('double_bond_radius', None)
    """ Double covalent-bond  radius. """
    self.triple_bond_radius       = kwargs.pop('triple_bond_radius', None)
    """ Triple covalent-bond  radius. """
    self.van_der_waals_radius       = kwargs.pop('van_der_waals_radius', None)
    """ van der Walls radius. """
    self.boiling_point       = kwargs.pop('boiling_point', None)
    self.melting_point       = kwargs.pop('melting_point', None)
    self.density             = kwargs.pop('density', None)
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
#   string = string.replace("alt\"", "alt=\"")
#   soup = BeautifulSoup(re.sub(re_swf,"rainbow.swf\"",string))

    atom = Element(name=name)
#   atom.symbol = soup.findChild( name="a", attrs={"title": "Element names and symbols"},\
#                                 text=" Symbol").parent.parent.contents[1].split()[1]
#   atom.atomic_number = soup.findChild(name="a", attrs={"title": "Element atomic numbers"})\
#                                      .parent.contents[-1].split()[1]
#   atom.atomic_number = int(atom.atomic_number)
#   atom.atomic_weight = soup.findChild(name="a", attrs={"title": "Element atomic weights"})\
#                                      .parent.prettify()
#   found = re_atomweight.search(atom.atomic_weight)
#   if found == None: print name
#   else: atom.atomic_weight = float(found.group(1))

    
#   # ionization stuff
#   if not exists(join("elements", name + "_atoms.html")):
#     file = urllib.urlopen("http://www.webelements.com/{0}/atoms.html".format(name.lower()))
#     string = file.read()
#     file.close()
#   else: 
#     with open(join("elements", name + "_atoms.html"), "r") as file: string = file.read()
#   soup = BeautifulSoup(string) 
#   # electron affinity
#   found = re.search("of\s+{0}\s+is\s+(\S+)".format(name.lower()), string)
#   if found.group(1) == "no": atom.electron_affinity = None
#   else: atom.electron_affinity = float(found.group(1)) * pq.kilo * pq.J / pq.mol
#   # ionization energies
#   energies = []
#   for child in soup.findChild(name="table", attrs={"class":"chemistry-data"})\
#                    .findChildren(name='td'):
#     energies.append(float(child.string) * pq.kilo * pq.J / pq.mol)
#   atom.ionization_energies = energies if len(energies) > 0 else None


#   # electronegativities.
#   if not exists(join("elements", name + "_electronegativity.html")):
#     file = urllib.urlopen("http://www.webelements.com/{0}/electronegativity.html"\
#                           .format(name.lower()))
#     string = file.read()
#     file.close()
#   else: 
#     with open(join("elements", name + "_electronegativity.html"), "r") as file:
#         string = file.read()
#   soup = BeautifulSoup(string) 
#   attrs = { "href": "../periodicity/electronegativity_pauling/",\
#             "title": "View definition and pictures showing periodicity "\
#                      "of Pauling electronegativity"}
#   pauling = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
#   pauling = pauling.split()[0]
#   atom.pauling = float(pauling) if pauling != "no" else None

#   attrs = { "href": "../periodicity/electronegativity_sanderson/" }
#   sanderson = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
#   sanderson = sanderson.split()[0]
#   atom.sanderson = float(sanderson) if sanderson != "no" else None

#   attrs = { "href": "../periodicity/electroneg_allred_rochow/" }
#   allred_rochow = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
#   allred_rochow = allred_rochow.split()[0]
#   atom.allred_rochow = float(allred_rochow) if allred_rochow != "no" else None

#   attrs = { "href": "../periodicity/electroneg_mulliken_jaffe/" }
#   mulliken_jaffe = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1]
#   if name in ["Germanium", "Gallium", "Carbon", "Lead", "Boron", "Silicon", "Tin",\
#               "Thallium", "Aluminium", "Indium"]: 
#     mulliken_jaffe = mulliken_jaffe.contents[0]
#   else: mulliken_jaffe = mulliken_jaffe.string
#   mulliken_jaffe = mulliken_jaffe.split()[0]
#   atom.mulliken_jaffe = float(mulliken_jaffe) if mulliken_jaffe != "no" else None

#   attrs = { "href": "../periodicity/electronegativity_allen/" }
#   allen = soup.findChild(name="a", attrs=attrs).parent.parent.contents[-1].string
#   allen = allen.split()[0]
#   atom.allen = float(allen) if allen != "no" else None
    
    # atom sizes
    if not exists(join("elements", name + "_atom_sizes.html")):
      file = urllib.urlopen("http://www.webelements.com/{0}/electronegativity.html"\
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
    atom.atomic_radius = float(atomic_radius) * pq.picometre if atomic_radius != "no" else None
    
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
      print name, atom.van_der_waals_radius
