
def single_site_factory(lattice, J0=False, J1=False, **mb):
  """ Creates clusters starting for a lattice with a single site. 

      Things are simpler for lattice with a single site. 
      Fortunately, one can transpose the results from a single site to those of a
      multisite.
  """ 
  from re import match
  from ..crystal import space_group, coordination_shells
  from ..error import ValueError
  from .cluster import Cluster, spin
  check_keywords(J0, J1, **mb)
  if len(lattice) != 1:
    raise ValueError("Expects lattice with a single site.")
  if not hasattr(lattice[0].type, '__len__'): 
    raise ValueError("Lattice occupation is singular.")
  if len(lattice[0].type) < 2: 
    raise ValueError("Lattice occupation is singular.")

  lattice = lattice.copy()
  lattice.space_group = space_group(lattice)
   
  result = []
  if J0: result.append(Cluster(lattice))
  if J1:
    result.append(Cluster(lattice))
    result[-1].spins = spin([0,0,0])
  # figure out max number of shells.
  maxshells = max(mb.keys(), key=lambda x: int(x[1:]))
  # compute these shells once and for all.
  shells = coordination_shells(lattice, maxshells, lattice[0].pos)

  # creates many bodies.
  for key in sorted(mb.keys()):
    if mb[key] <= 0: continue
    order = int(key[1:])
    dummy = _create_clusters(lattice, shells=shells[:mb[key]], order=order)
    result.extend(dummy)
  return result

def check_keywords(J0, J1, **mb):
  # checks that keywords are well formed.
  from re import match, compile
  from ..error import KeyError, ValueError
  key_regex = compile("B(\d+)")
  for key in mb.keys(): 
    a_re = match(key_regex, key)
    if a_re is None:
      raise KeyError("Keyword {0} is not of the form B(\d+)".format(key))
    if a_re.end() != len(key):
      raise KeyError("Keyword {0} is not of the form B(\d+)".format(key))
    if int(a_re.group(1)) < 2:
      raise KeyError("Bn with n < 2.")
    if mb[key] < 0:
      raise ValueError("Negative coordination shell for {0}.".format(key))
  if len(mb) == 0 and (not J0) and (not J1): 
    raise ValueError("No clusters to create.")

def _create_clusters(lattice, shells, order):
  """ Creates a set of clusters of given order up to given shell. """
  from numpy import array, any
  from .cluster import Cluster
  from .cppwrappers import ProductILJ
  # create lists of neighbors from shells.
  neighbors = array([u[1] for v in shells for u in v])
  # order them same as clusters.
  knowns = set()
  results = []
  # loop over figures made from those neighbors.
  # each figure is represented by an ordered set of indices in the neighbor
  # list.  This makes it easy to recognize which figures are already present in
  # the list of clusters.
  for indices in ProductILJ(range(len(neighbors)), repeat=order-1):
    # make sure this figure was not previously found by symmetry.
    if indices in knowns: continue
    # now create cluster
    a = Cluster()
    a.spins = array( [ (neighbors[i], 0) for i in indices], 
                     dtype=[('position', 'f8', 3), ('sublattice', 'i8', 1)] )
    a.add_spin([0,0,0])
    results.append(a)
    
    # now we add to the list of known figures those which are symmetrically
    # equivalent to this one.
    a._create_symmetrized()
    for u in a._symmetrized: 
      indices = []
      for spin in a.spins['position'][1:]:
        for i, neigh in enumerate(neighbors):
          if not any(abs(spin - neigh) > 1e-8):
            indices.append(i)
            break;
      indices = sorted(indices)
      if order == 2 or all(u < v for u, v in zip(indices[:-1], indices[1:])):
        knowns.add(set(indices))
  return results
