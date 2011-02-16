""" Package to enumerate structural decorations of a lattice. """
__docformat__ = "restructuredtext en"
__all__ = ['Enum', 'as_structure', 'as_numpy']
from _enumeration import as_structure, as_numpy   
from ..crystal import Lattice

def self_as_global_lattice(method):
  """ Sets the global lattice to self on entering function. """
  def wrapper(self, *args, **kwargs):
    from ..crystal import Structure
    oldlattice = None
    try: oldlattice = Structure().lattice
    except: pass
    self.set_as_crystal_lattice()
    result = method(self, *args, **kwargs)
    if oldlattice != None: oldlattice.set_as_crystal_lattice()
    return result
  wrapper.__name__   = method.__name__
  wrapper.__doc__    = method.__doc__
  wrapper.__module__ = method.__module__
  wrapper.__dict___  = method.__dict__.update(method.__dict__)
  return wrapper

    
class Enum(Lattice):
  """ Enhanced `crystal.Lattice` providing generators over inequivalent configurations. """

  def __init__(self, lattice=None, tolerance=1e-12):
    """ Initializes the enumeration class. 
    
        :Parameters: 
          lattice : `crystal.Lattice`
            If this parameter is provided, then the lattice parameters are deep
            copied.
          tolerance : float
            Tolerance with which to find crystal symmetries. 
    """
    Lattice.__init__(self)
    if lattice == None: return

    from copy import deepcopy
    # Copies information from lattice if provided.
    self.cell = lattice.cell.copy()
    self.scale = lattice.scale
    self.name = lattice.name
    for site in lattice.sites: 
      self.sites.append(deepcopy(site))
    self.__dict__.update(deepcopy(lattice.__dict__))
    self.find_space_group()

  @property 
  def nflavors(self):
    """ Number of flavors in this lattice. """
    from _enumeration import count_flavors
    return count_flavors(self)

  @property
  def nsites(self):
    """ Number of sites in lattice. """
    return len([0 for u in self.sites if  len(u.type) > 1])

  @property
  def transforms(self):
    """ Transforms for this particular lattice. """
    from _enumeration import create_transforms
    return create_transforms(self)

  def supercells(self, n):
    """ Iterates over supercells. """
    from _enumeration import find_all_cells
    for cell in find_all_cells(self, n): yield cell

  def smiths(self, n):
    """ Iterates over smith groups. """
    from _enumeration import find_all_cells, create_smith_groups
    supercells = find_all_cells(self, n)
    for smith in create_smith_groups(self, supercells): yield smith

  @self_as_global_lattice
  def xn(self, n, callback = None):
    """ Iterates over all decorations with n atoms. 
    
        :Parameters:
          n : int
            number of atoms to consider.
          callback : callable or None
            Callable taking x, flavorbase, and smith on call. 
            It should return a tuple where the first item says to proceed with
            the loop, and the other wether to set that entry in the database to
            True or False. If the second item is None, then the database entry
            is untouched.

         :return: Yields a tuple consisting of x, the smith structure,
                  the supercell, and the flavorbase.
    """
    from numpy import dot
    from _enumeration import LabelExchange, create_flavorbase, Translation, Database, as_numpy
    assert len(self.sites) > 0, ValueError("No sites in lattice.")
    nsites = self.nsites
    assert nsites > 0, ValueError("Current lattice does not have sites with multiple occupation.")
    nflavors = self.nflavors
    transforms = self.transforms

    # loop over smith groups.
    for smith in self.smiths(n):
      card           = int(smith.smith[0]*smith.smith[1]*smith.smith[2]*nsites)
      label_exchange = LabelExchange(card, nflavors )
      flavorbase     = create_flavorbase(card, nflavors)
      translations   = Translation(smith.smith, nsites)
      database       = Database(card, nflavors)
      maxterm        = 0

      # loop over possible decorations.
      for x in xrange(1, int(pow(nflavors, card))-1):
 
        # loops now if false.
        if not database[x]: continue

        # callback
        if callback != None:
          loop, entry = callback(x, flavorbase, smith)
          if entry != None: database[x] = entry
          if loop: continue

        # check for supercell independent transforms.
        maxterm = x
        # loop over permutations
        for labelperm in label_exchange:
          t = labelperm(x, flavorbase)
          if t > x: database[t] = False
        # loop over translational symmetries.
        for translation in translations:
          t = translation(x, flavorbase)
          if t > x: database[t] = False
          elif t == x: 
            database[t] = False
            continue
          for labelperm in label_exchange:
            u = labelperm(t, flavorbase)
            if u > x: database[u] = False

    # checks supercell dependent transforms.
    for nsupercell, supercell in enumerate(smith.supercells):
      mine = []
      # creates list of transformation which leave the supercell invariant.
      cell = dot(self.cell, supercell.hermite)
      specialized = []
      for transform in transforms:
        if not transform.invariant(cell): continue
        transform.init(supercell.transform, smith.smith)
        if not transform.is_trivial: specialized.append( transform )

      specialized_database = Database(database)
      for x in xrange(1, maxterm+1):
        if not database[x]: continue
        maxterm = x
        
        for transform in specialized:
          t = transform(x, flavorbase)
          if t == x: continue
          specialized_database[t] = False

          for labelperm in label_exchange:
            u = labelperm(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            
          for translation in translations:
            u = translation(t, flavorbase)
            if u == x: continue
            specialized_database[u] = False
            for labelperm in label_exchange:
              v = labelperm(u, flavorbase)
              if v == x: continue
              specialized_database[v] = False
        if specialized_database[x]: yield x, smith, supercell, flavorbase

  def structures(self, *args, **kwargs):
    """ Yields inequivalent structures.
    
        Parameters are those of `Enum.xn`.
        This generator yields actual `crystal.Structure.`
    """
    from numpy import zeros
    from _enumeration import as_structure
    oldhermite = zeros((3,3))
    for x, smith, supercell, flavorbase in self.xn(*args, **kwargs):
      hermite = supercell.hermite
      if oldhermite != hermite: structure = self.to_structure(dot(self.cell, hermite))
      yield as_structure(structure, x, flavorbase)
