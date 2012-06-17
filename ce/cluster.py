def spin(position=None, flavor=0, sublattice=0, subgroup=0):
  """ Returns a spin. 

      A spin is a numpy `structured record`_  with fields for the position,
      flavor, and sublattice.
  """ 
  from numpy import array
  if position is None: position = 0,0,0
  return array( [position, flavor, sublattice], 
                dtype=[ ('position', 'f8', 3), 
                        ('flavor', 'i8', 1), 
                        ('sublattice', 'i8', 1),
                        ('subgroup', 'i8', 1) ] )

class Cluster(object):
  def __init__(self, lattice):
    from ..crystal import Structure, spacegroup
    from ..error import ValueError

    super(Cluster, self).__init__()
    if not isinstance(lattice, Structure): 
      raise ValueError( "Lattice should be an instance of "                    \
                        "lada.crystal.Structure.")
    self.lattice = lattice
    """ Lattice back-bone on which this cluster exists. 

        Should be an instance of
        :py:class:`~lada.crystal.cppwrappers.Structure` or derived.
    """
    self.spins = None
    """ Holds array of spins. """
    self.spacegroup = spacegroup(self.lattice)
    """ Symmetry operations. """
    self._sitesmap = self._count_sites()
    """ Map from CE sublattice to lattice sites. """


  @property
  def nbsublattices(self): 
    """ Number of sublattices. """
    return len(self._sitesmaps)
  @property
  def order(self):
    """ Returns order of figure, eg number of sites involved. """
    return 0 if self.spins is None else len(self.spins)
  def occupation(self, subgroup=0):
    """ Occupation of a particular subgroup of atoms.
    
        A subgroup contains all equivalent atoms in a multicomponent lattice.
    """
    from ..error import IndexError
    if subgroup <= -len(self._sitesmap) or subgroup >= len(self._sitesmap): 
      raise IndexError('Incorrect subgroup index.')
    return set(self._sitesmap[subgroup][0].type)
    
  def add_spin(self, position, flavor=0):
    """ Adds a spin to array. 
    
        Changes array of spins. Hence, one should not keep a reference to a
        cluster's array of spins until it is complete.
    """
    from numpy import concatenate
    from ..crystal import which_site
    from ..error import IndexError, ValueError
    # remove symmetrization.
    if hasattr(self, '_symmetrized'): del self._symmetrized

    # find sublattice
    sublattice = which_site(position, self.lattice)
    if sublattice == -1: raise ValueError('Could not find atomic site')
    site = self.lattice[sublattice]

    # finds subgroup
    found = False
    for subgroup, group in enumerate(self._sitesmap):
      if id(site) in [id(u) for u in group]: 
        found = True
        break
    if found == False: 
      raise ValueError('Sublattice position is not part of CE.')
      
    if flavor <= len(site.type) or flavor >= len(site.type):
      raise IndexError('Incorrect flavor for given sublattice.')

    # create spin and adds it to group.
    spins = [spin(position, flavor, sublattice, subgroup)]
    self.spins = spins if self.spins is None                                   \
                 else concatenate((self.spins, spins), axis=0)

  def _get_symmetrized(self):
    """ Computes equivalent clusters. """
    self._symmetrized = []
    for op in self.spacegroup:
      # transform spins
      cluster = self.__class__()
      for spin in self.spins:
        cluster.add_spin(op[:3] * spin['position'] + op[3], spin['flavor'])
      # check if exists in symmetrized bits.
      if cluster in self._symmetrized: continue 


  def _analyze_lattice(self): 
    """ Analyzes lattice. 

        Can only be called once. 
    """
    from ..error import internal
    from ..crystal import equivalence_iterator
    if hasattr(self, '_nbsites'): 
      raise internal('_analyze_lattice cannot be called twice.')
    def split_occ(a, b):
      """ split along occupations. """
      hasiter_a = hasattr(a.type, '__iter__')
      hasiter_b = hasattr(b.type, '__iter__')
      if (hasiter_a and hasiter_b): 
        return set(a.type) == set(b.type)
      if (not hasiter_a) and (not hasiter_b): return a == b
      return False
    equivsg = equivalence_iterator( self.lattice,                              \
                                    operators=self.spacegroup,                 \
                                    splitocc=split_occ )
    equivsg = [list(u) for u in equivsg]
    return [u for u in equivsg if len(u) > 1]
