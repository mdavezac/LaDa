def cluster_factory(lattice, J0=False, J1=False, **mb):
  """ Returns class of cluster classes. 

      Creates a list of clusters using a simple input. 

      :param lattice:
          The backbone lattice on which the Ising model is created.
      :param bool J0:
          If True, adds zero order term.
      :param bool J1: 
          If True, adds on-site terms.
      :param mb: 
          All other many body terms, where the keys should be "B2", "B3", "Bn",
          where n is the order of the many-body interaction, and the values is
          the number of shells to look for.
  """
  from numpy import dot
  from ._single_site_factory import check_keywords, factory as ss_factory
  from .cluster import Cluster, spin
  from ..crystal import which_site, space_group, coordination_shells

  # checks that keywords are well formed.
  check_keywords(J0, J1, **mb)

  if len(lattice) == 1: return ss_factory(lattice, J0, J1, **mb)

  # add space group if it does not already exist.
  if not hasattr(lattice, 'spacegroup'): 
    lattice.spacegroup = space_group(lattice)

  # computes equivalent sites and add stuff to lattice for later
  nonequiv = set( i for i in range(len(lattice)) ) 
  for i, site in enumerate(lattice):
    if i not in nonequiv: continue
    for op in lattice.spacegroup:
      j = which_site( dot(op[:3], site.pos) + op[3], lattice )
      if j != i and j in nonequiv: nonequiv.remove(j)
  lattice = lattice.copy()
  for i, site in enumerate(lattice):
    site.asymmetric = i in nonequiv
    site.index = i

  result = []
  # First creates J0
  if J0:
    result.append(Cluster(lattice))
  # then creates J1 for these sites.
  if J1: 
    for i in nonequiv: 
      site = lattice[i]
      if not hasattr(site.type, '__iter__'): continue
      result.append(Cluster(lattice))
      result[-1].spins = spin([0,0,0], site.index)

  if len(mb) == 0: return result
  # figure out max number of shells.
  maxshells = max(mb.itervalues())
  # creates many bodies.
  for i in nonequiv:
    site = lattice[i]
    if not hasattr(site.type, '__iter__'): continue
    # creates list of neighbors with indices into shells.
    allshells = coordination_shells(lattice, maxshells, site.pos)
    shells = []
    shellindices = []
    for shell in allshells:
      shellindices.append(len(shell))
      for atom, vector, distance in shell:
        if not hasattr(atom.type, '__iter__'): continue
        if not hasattr(atom.type, '__len__'): continue
        if len(atom.type) < 2: continue
        pos = vector.copy()
        if atom is not site: pos += site.pos - atom.pos
        shells.append(spin(pos, atom.index))
    # loop over manybody figures and creates them.
    for key in sorted(mb.keys()):
      if mb[key] <= 0: continue
      order = int(key[1:])
      dummy = _create_clusters( lattice, 
                                neighbors=shells[:shellindices[mb[key]-1]],
                                site=site,
                                order=order )
      result.extend(dummy)
  return result

def _create_clusters(lattice, neighbors, site, order):
  """ Creates a set of clusters of given order up to given shell. """
  from numpy import any, concatenate, not_equal
  from .cluster import Cluster, spin
  from .cppwrappers import ProductILJ
  knowns = set()
  results = []
  # loop over figures made from those neighbors.
  # each figure is represented by an ordered set of indices in the neighbor
  # list.  This makes it easy to recognize which figures are already present in
  # the list of clusters.
  for indices in ProductILJ(range(len(neighbors)), order-1):
    # make sure this figure was not previously found by symmetry.
    if indices in knowns: continue
    # now create cluster
    a = Cluster(lattice)
    a.spins = concatenate( [spin([0,0,0], site.index)] 
                           + [neighbors[i] for i in indices] )
    results.append(a)
    
    # now we add to the list of known figures those which are symmetrically
    # equivalent to this one.
    a._create_symmetrized()
    for u in a._symmetrized: 
      if u[0]['sublattice'] != site.index: continue
      indices = []
      for vector in u[1:]:
        for i, neigh in enumerate(neighbors):
          if not any(not_equal(vector, neigh[0])):
            indices.append(i)
            break;
      indices = tuple(sorted(indices))
      if order == 2 or all(u < v for u, v in zip(indices[:-1], indices[1:])):
        knowns.add(indices)
  return results
