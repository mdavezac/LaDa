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
  from re import compile, match
  from .cluster import Cluster
  from ..crystal import which_site, space_group

  # checks that keywords are well formed.
  key_regex = compile("B(\d+)")
  for key in mb.keys(): 
    a_re = match(key_regex, key)
    assert a_re is not None, "Keyword %s is not of the form B(\d+)" % (key)
    assert a_re.end() == len(key), "Keyword %s is not of the form B(\d+)" % (key)
    assert int(a_re.group(1)) > 1, "Cannot understand keyword %s" % (key)
    assert mb[key] >= 0, "Cannot understand input %s=%i" % (key, mb[key])

  # sanity check.
  assert len(mb) > 0 or J0 or J1, ValueError("No clusters to create.")
  lattice.space_group = space_group(lattice)

  # computes equivalent sites and add stuff to lattice for later
  nonequiv = set( i for i in range(len(lattice)) ) 
  for i, site in enumerate(lattice):
    if i not in nonequiv: continue
    for op in lattice.space_group:
      j = which_site( dot(op[:3], site.pos) + op[3], lattice )
      if j != i and j in nonequiv: nonequiv.remove(j)
  lattice = lattice.copy()
  for i, site in enumerate(lattice):
    site.asymmetric = i in nonequiv
    site.index = i

  result = []
  # now creates multi-lattice clusters, along with index bookkeeping.
  # First creates J0
  if J0:
    result.append(Cluster(lattice))
  # then creates J1 for these sites.
  if J1: 
    for i in nonequiv: 
      site = lattice[i]
      if not hasattr(site.type, '__iter__'): continue
      for flavor in xrange(len(site.type)-1): 
        cluster = Cluster(lattice)
        cluster.add_spin(site.pos, flavor)
        result.append(cluster)
  # creates many bodies.
  for key in sorted(mb.keys()):
    if mb[key] <= 0: continue
    regex = match(key_regex, key)
    order = int(regex.group(1))
    dummy = _create_clusters(lattice, nshells=mb[key], order=order)
    result.extend(dummy)
  return result

def _create_clusters(lattice, nshells, order):
  """ Creates sets of clusters of given order up to given shell """
  from itertools import product, combinations, chain, imap,                     \
                        combinations_with_replacement as cwr
  from copy import deepcopy
  from numpy import concatenate
  from ..crystal import coordination_shells
  from .cluster import Cluster, spin

  def pick_loop(which_shells, shells, site):
    """ Flattens loop over positions and flavors """
    ws = set(which_shells)
    iterables = [ (xrange(len(shells[s])), which_shells.count(s))              \
                  for s in set(which_shells) ]
    iterables = [combinations(l, c) for l, c in iterables]
    for indices in product(*iterables):
      indices = _flatten(*indices)
      positions = [ shells[which_shells[i]][index] 
                    for i, index in enumerate(indices) ]
      yield ( ( p[1] if p[0] is site else p[1] - p[0].pos + site.pos,          \
                p[0].index )                                                   \
              for p in positions )
      continue
      
  firsts = [atom for atom in lattice if atom.asymmetric]
  figures = []
  for sublattice, site in enumerate(firsts): 
    if not hasattr(site.type, '__iter__'): continue
    if len(site.type) < 2: continue
    # compute coordination shells, removing those atoms which do not have the
    # more than one type
    shells = []
    for shell in coordination_shells(lattice, nshells, site.pos):
      intermediates = [p for p in shell if hasattr(p[0].type, '__iter__')]
      intermediates = [p for p in intermediates if len(p[0].type) > 1]
      if len(intermediates) > 0: shells.append(intermediates)
    
    # Loop over combination of shells. 
    # We pick spins explicitely from the coordination shells.
    # Since each shell contains all symmetrically equivalent spins (possibly
    # more than one group of symmetrically equivalent spins, but this is
    # easier), we need to check symmetry issues within the combination
    # loop preferentially (eg the one right below this statement). 
    for whichshells in cwr(range(len(shells)), order-1):
      intermediates = []
      # check we are not asking for too many spins from same shell.
      if any( len(shells[s]) < whichshells.count(s) for s in whichshells ): 
        continue
      # This loops over blocks of atomic positions and associated flavor.
      # Within this loop, the clusters can be created and tested.
      for block in pick_loop(whichshells, shells, site):
        block = [k for k in block]
        print site.index, block
        a = Cluster(lattice)
        iterable = chain([([0,0,0], sublattice)], block) 
        a.spins  = concatenate([spin(*u) for u in iterable])
        if any(b._contains(a.spins) for b in intermediates): continue
        if any(b._contains(a.spins) for b in figures): continue
        a._create_symmetrized()
        intermediates.append(a)
      figures.extend(intermediates)

  results = []
  for r in figures:
    spins = r.spins
    iterable = (range(len(lattice[s].type)-1) for s in spins['sublattice'])
    intermediates = []
    for flavors in product(*[k for k in iterable]):
      spins['flavor'] = flavors
      if any(b._contains(spins) for b in intermediates): continue
      a = Cluster(lattice, spins.copy())
      a._create_symmetrized()
      intermediates.append(a)
    results.extend(intermediates)
  return figures


def _flatten(*args):
  """ Flatten arguments """
  result = []
  for arg in args:
    if hasattr(arg, '__iter__'): result.extend(_flatten(*arg))
    else:                        result.append(arg)
  return result

 



