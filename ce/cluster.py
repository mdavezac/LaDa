def spin(position=None, sublattice=0, flavor=0):
  """ Returns a spin. 

      A spin is a numpy `structured record`_  with fields for the position,
      flavor, and sublattice.
  """ 
  from numpy import array
  if position is None: position = 0,0,0
  return array( [(position, sublattice, flavor)], 
                dtype=[ ('position', 'f8', 3), 
                        ('sublattice', 'i8', 1), 
                        ('flavor', 'i8', 1) ] )

class Cluster(object):
  def __init__(self, lattice, spins=None):
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
    self.spins = spins
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

    # find sublattice
    sublattice = which_site(position, self.lattice)
    if sublattice == -1: raise ValueError('Could not find atomic site')
    site = self.lattice[sublattice]

    # keep negative flavors as markers
    if flavor >= 0 and flavor >= len(site.type) - 1:
      raise IndexError( 'Incorrect flavor for given sublattice.\n'            \
                        '   0 <= flavor < {0}\n'.format(len(site.type)-1) )

    # create spin and adds it to group.
    return spin(position - self.lattice[sublattice].pos, sublattice, flavor)

  def add_spin(self, position, flavor=0):
    """ Adds a spin to array. 
    
        Changes array of spins. Hence, one should not keep a reference to a
        cluster's array of spins until it is complete.

        :param position:
          Atomic site on the back-bone lattice.
        :type position: 3d-vector
        :param (int) flavor:
          Flavor of the function for the particular atomic-site. It should be
          :math:`0 \leq \text{flavor} < M-1`, where :math:`M` is the number of
          possible species on the atomic site. This defaults to 0 and can be
          ignored in the case of binary cluster expansions.
    """
    from operator import itemgetter
    from numpy import concatenate, any, all, sum, logical_and
    from ..error import ValueError
    # remove symmetrization.
    if hasattr(self, '_symmetrized'): del self._symmetrized
    # create spin
    spin = self._spin(position, flavor)
    # if first spin, easy way out.
    if self.spins is None:
      self.spins = spin
      return 
    # otherwise, check whether it is already in the cluster
    positions   = all(self.spins['position'] == spin['position'], axis=1) 
    sublattices = self.spins['sublattice'] == spin['sublattice']
    if any(logical_and(positions, sublattices)): 
      raise ValueError('Spin already present in cluster.')
    # Sort the spins. Since we also require that at least one spin is in the
    # unit-cell, this should mean fewer calculations at the  end of the day.
    self.spins = concatenate((self.spins, spin), axis=0)
    spins = sum(self.spins['position']* self.spins['position'], axis=1)
    spins = [ (s, l, f, p[0], p[1], p[2])                                      \
              for s, l, f, p in zip( spins, self.spins['sublattice'], 
                                     self.spins['flavor'],
                                     self.spins['position']) ]
    spins = [u[0] for u in sorted(enumerate(spins), key=itemgetter(1))]
    self.spins = self.spins[spins].copy()
    
  def _contains(self, spins):
    """ Compares cluster to those already in _symmetrized """
    from numpy import all, any
    if len(spins) != self.order: return False
    for cluster in self._symmetrized: 
      if all([any([all(v == s) for s in spins]) for v in cluster]):
        return True
    return False

  def _create_symmetrized(self):
    """ Computes equivalent clusters. """
    from numpy import concatenate, all, dot, abs
    from ..crystal import into_voronoi as iv
    from ..error import ValueError
    if self.spins is None or len(self.spins) == 0: return
    if all(abs(self.spins['position'][0]) > 1e-8):
      raise ValueError('Expected first spin in the unit cell.')
    self._symmetrized = []
    positions = [ p + self.lattice[s].pos                                      \
                  for p, s in zip( self.spins['position'], 
                                   self.spins['sublattice'] ) ]
    for op in self.spacegroup:
      # transform spins
      spins = [ self._spin(dot(op[:3],p)+op[3], flavor)                        \
                for p, flavor in zip(positions, self.spins['flavor']) ]
      spins = concatenate(spins)
      assert all(abs(iv(spins[0]['position'], self.lattice.cell)) < 1e8)
      spins['position'] -= spins[0]['position']
      # check if exists in symmetrized bits.
      if not self._contains(spins): self._symmetrized.append(spins)

  def occupation_mapping(self):
    """ Map of occupations to cluster function
    
        A mapping which returns the value of a particular spin according to
        occupation and flavor:
          
        .. code-block::

          mapping = cluster.occupation_mapping()
          mapping[atom.site][spin_flavor][atom.type]

        where ``atom`` is an element of a structure, ``atom.site`` is likely
        the index of the atomic site into the backbone lattice, and atom.
        ``flavor`` is the flavor of a spin this particular cluster.

        This particular mapping is taken from Ref.[*]_.

        .. [*]: `Multicomponent multisublattice alloys, nonconfigurational
                 entropy and other additions to the Alloy Theoretic Automated
                 Toolkit`__, Axel van de Walle, Calpha, *33*, 266-278 (2009)
        .. __: http://dx.doi.org/10.1016/j.calphad.2008.12.005
    """
    from numpy import cos, sin, pi
    result = []
    for site in self.lattice:
      if not hasattr(site.type, '__iter__'): result.append(None); continue
      if len(site.type) < 2: result.append(None); continue
      result.append([])
      for flavor in xrange(len(site.type)-1): 
        result[-1].append({})
        for i, type in enumerate(sorted(site.type)):
          arg = 2e0 * pi*float(flavor+1) * int(float(i) * 0.5 + 0.5)           \
                / float(len(site.type))
          if flavor % 2 == 0: result[-1][-1][type] = -cos(arg)
          else: result[-1][-1][type] = -sin(arg)
    return result

  def __call__(self, structure, mapping=None, transform=None, fhmapping=None):
    """ Returns PI of a particular structure. """
    from lada.crystal import HFTransform
    if self.spins is None or len(self.spins) == 0:
      return float(len(structure))
    if not hasattr(self, '_symmetrized'): self._create_symmetrized()
    if mapping is None: mapping = self.occupation_mapping()
    if transform is None: transform = HFTransform(self.lattice, structure)
    if fhmapping is None: fhmapping = fhmap(self.lattice, structure)
    result = 0e0
    # loop over possible origin of cluster
    for origin in structure:
      pos = origin.pos - self.lattice[origin.site].pos
      # loop over symmetrically equivalent clusters
      for spins in self._symmetrized:
        # no need to double count
        if spins['sublattice'][0] != origin.site: continue
        intermediate = 1e0
        # loop over each spin in cluster
        for spin in spins: 
          index = transform.index(pos + spin['position'], spin['sublattice'])
          atom  = structure[ fhmapping[index] ]
          intermediate *= mapping[atom.site][spin['flavor']][atom.type]
        result += intermediate
    return result

  def __str__(self):
    if self.order == 0: return 'J0'
    if self.order == 1:
      sublat = self.spins['sublattice'][0]
      flavor = self.spins['flavor'][0]
      nflavr = -1 if flavor < 0 else len(self.lattice[sublat].type) - 1
      return 'J1: site {0}, flavor {1}/{2}'.format(sublat, flavor, nflavr)
    result = 'M{0}:\n{1}\n'.format(self.order, self.spins)
    return result


def fhmap(lattice, supercell):
  """ Creates a map from Forkade-Hart indices to atom indices """
  from lada.error import ValueError
  from lada.crystal import HFTransform
  transform = HFTransform(lattice, supercell)
  results = {}
  for i, atom in enumerate(supercell):
    if not hasattr(atom, 'site') and len(lattice) > 0: 
      raise ValueError( 'Atoms in supercell were not mapped to lattice. '      \
                        'Please use lada.crystal.map_sites.' )
    site = getattr(atom, 'site', 0)
    try: index = transform.index(atom.pos - lattice[site].pos, site)
    except ValueError:
      print supercell.cell
      print atom.pos, lattice[site].pos, atom.pos-lattice[site].pos
      raise
    results[index] = i
  return results
