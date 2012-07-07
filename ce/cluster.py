def spin(position=None, flavor=0, sublattice=0):
  """ Returns a spin. 

      A spin is a numpy `structured record`_  with fields for the position,
      flavor, and sublattice.
  """ 
  from numpy import array
  if position is None: position = 0,0,0
  return array( [(position, flavor, sublattice)], 
                dtype=[ ('position', 'f8', 3), 
                        ('flavor', 'i8', 1), 
                        ('sublattice', 'i8', 1) ] )

class Cluster(object):
  def __init__(self, lattice):
    from lada.crystal import Structure, space_group
    from lada.error import ValueError

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
    self.spacegroup = space_group(self.lattice)
    """ Symmetry operations. """

  @property
  def nbsublattices(self): 
    """ Number of sublattices. """
    return len(self._sitesmaps)
  @property
  def order(self):
    """ Returns order of figure, eg number of sites involved. """
    return 0 if self.spins is None else len(self.spins)
    
  def _spin(self, position, flavor=0):
    """ Returns formatted spin. """
    from lada.crystal import which_site
    from lada.error import IndexError, ValueError
    # remove symmetrization.
    if hasattr(self, '_symmetrized'): del self._symmetrized

    # find sublattice
    sublattice = which_site(position, self.lattice)
    if sublattice == -1: raise ValueError('Could not find atomic site')
    site = self.lattice[sublattice]

    if flavor <= len(site.type) or flavor >= len(site.type):
      raise IndexError('Incorrect flavor for given sublattice.')

    # create spin and adds it to group.
    return spin(position, flavor, sublattice)

  def add_spin(self, position, flavor=0):
    """ Adds a spin to array. 
    
        Changes array of spins. Hence, one should not keep a reference to a
        cluster's array of spins until it is complete.
    """
    from numpy import concatenate
    spin = self._spin(position, flavor)
    if self.spins is None: self.spins = spin
    elif any(self.spins == spin): 
      raise ValueError('Spin already present in cluster.')
    else: self.spins = concatenate((self.spins, spin), axis=0)
    
  def _contains(self, spins):
    """ Compares cluster to those already in _symmetrized """
    from numpy import all, any
    if len(spins) != self.order: return False
    for cluster in self._symmetrized: 
      if all(any(cluster.spins == s) for s in spins): return True
    return False


  def _get_symmetrized(self):
    """ Computes equivalent clusters. """
    from numpy import concatenate
    self._symmetrized = []
    for op in self.spacegroup:
      # transform spins
      spins = [ self._spin(op[:3] * spin['position'] + op[3], spin['flavor'])
                for spin in self.spins ]
      spins = concatenate(spins)
      # check if exists in symmetrized bits.
      if not self._contains(spins): self._symmetrized.append(spins)
