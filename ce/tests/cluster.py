def get_cell(n=5):
  from numpy.random import randint
  from numpy.linalg import det
  cell = randint(2*n, size=(3,3)) - n
  while abs(det(cell)) < 1e-8:
    cell = randint(2*n, size=(3,3)) - n
  return cell

def test_addspins():
  """ Test adding spins to cluster. 
     
      Check failure modes.
  """ 
  from numpy import all
  from lada.crystal import binary
  from lada.ce import Cluster
  from lada.ce.cluster import spin
  from lada.error import ValueError, IndexError

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = Cluster(lattice)
  
  # Wrong position
  try: a.add_spin(lattice[0].pos+0.1)
  except ValueError: pass
  else: raise Exception()
  # Wrong flavor
  try: a.add_spin(lattice[0].pos, 1)
  except IndexError: pass
  else: raise Exception()

  # now add first spin
  a.add_spin(lattice[0].pos)
  assert len(a.spins) == 1
  assert all(a.spins[0] == spin([0, 0, 0], 0, 0))

  # try adding it again
  try: a.add_spin(lattice[0].pos, 0)
  except ValueError: pass
  else: raise Exception()

  # Wrong position
  try: a.add_spin(lattice[1].pos+[1.1, -0.5, -2.5])
  except ValueError: pass
  else: raise Exception()
  # Wrong flavor
  try: a.add_spin(lattice[0].pos+[1.0, -0.5, -2.5], 2)
  except IndexError: pass
  else: raise Exception()

  # Then add a differrent spin
  a.add_spin(lattice[1].pos + [1.0, -0.5, -2.5], 1)
  assert len(a.spins) == 2
  assert all(a.spins[0] == spin([0, 0, 0], 0, 0))
  assert all(a.spins[1] == spin([1.0, -0.5, -2.5], 1, 1))
  # position already added, different flavor
  try: a.add_spin(lattice[1].pos + [1.0, -0.5, -2.5], 0)
  except ValueError: pass
  else: raise Exception()

def test_spins_are_sorted():
  """ Check spin sorting when using add_spins. """
  from numpy import all, dot
  from numpy.random import randint
  from random import choice
  from itertools import permutations
  from lada.ce.cluster import spin
  from lada.ce import Cluster
  from lada.crystal import binary
  
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  # Trial with known result.
  a = Cluster(lattice)
  a.add_spin(lattice[1].pos, 1)
  a.add_spin(lattice[0].pos + [0.5, 0, 0.5], 0)
  a.add_spin(lattice[1].pos + [-0.5, 0, 0.5], 0)
  a.add_spin(lattice[0].pos + [-2.5, -1.0, 0.5], 0)
  assert all(a.spins[0] == spin([0, 0, 0], 1, 1))
  assert all(a.spins[1] == spin([0.5, 0, 0.5]))
  assert all(a.spins[2] == spin([-0.5, 0, 0.5], 1, 0))
  assert all(a.spins[3] == spin([-2.5, -1.0, 0.5], 0, 0))
  spins = a.spins.copy()
  for p in permutations(range(len(spins))):
    a = Cluster(lattice)
    for i in p:
      a.add_spin( spins[i]['position'] + lattice[spins[i]['sublattice']].pos,
                  spins[i]['flavor'] )
    assert all(a.spins == spins)

  # Trial with unknown result.
  for i in xrange(20):
    a = Cluster(lattice)
    for j in xrange(5): 
      site = choice([0, 1])
      flavor = 0
      if site == 1: flavor = choice([0, 1])
      pos = lattice[site].pos + dot(lattice.cell, randint(4, size=(3,))-2)
      try: a.add_spin(pos, flavor)
      except ValueError: pass
    if a.order < 1: continue
    spins = a.spins.copy()
    for p in permutations(range(len(spins))):
      a = Cluster(lattice)
      for i in p:
        a.add_spin( spins[i]['position'] + lattice[spins[i]['sublattice']].pos,
                    spins[i]['flavor'] )
      assert all(a.spins == spins)

def test_cmp():
  """ Test Cluster._contains function """
  from numpy import all, any
  from lada.crystal import binary
  from lada.ce import Cluster
  from lada.ce.cluster import spin

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  def cmp(a,b):
    if len(a) != len(b): return False
    return all([any([all(v == s) for s in a]) for v in b])

  a = Cluster(lattice)
  a.add_spin(lattice[0].pos)
  assert cmp(a.spins,     [spin([0, 0, 0], 0, 0)])
  assert not cmp(a.spins, [spin([0.5, 0.5, 0.5], 0, 0)])
  assert not cmp(a.spins, [spin([0, 0, 0], 0, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0], 1, 0)])
  assert not cmp(a.spins, [spin([0, 0, 0], 1, 1)])
  a = Cluster(lattice)
  a.add_spin(lattice[1].pos, 1)
  assert cmp(a.spins,     [spin([0, 0, 0], 1, 1)])
  assert not cmp(a.spins, [spin([0.5, 0.5, 0.5], 1, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0], 0, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0], 1, 0)])
  assert not cmp(a.spins, [spin([0, 0, 0], 0, 0)])
  a.add_spin(lattice[0].pos)
  assert cmp(a.spins,     [spin([0, 0, 0], 1, 1), spin([0, 0, 0])])
  assert cmp(a.spins,     [spin([0, 0, 0]), spin([0, 0, 0], 1, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0], 1, 1), spin([0, 0, 0]), spin([1, 0, 0])])
  assert not cmp(a.spins, [spin([1, 0, 0]), spin([0, 0, 0], 1, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0], 1), spin([0, 0, 0], 1, 1)])
  assert not cmp(a.spins, [spin([0, 0, 0]), spin([0, 0, 0], 1, 0)])
  
def test_occmap():
  from numpy import cos, sin, pi
  from lada.crystal import binary
  from lada.ce import Cluster

  lattice = binary.zinc_blende()
  for atom in lattice: atom.type = ['Si', 'Ge']

  a = Cluster(lattice)
  mapping = a.occupation_mapping()
  assert len(mapping) == len(lattice)
  assert len(mapping[0]) == 1
  assert len(mapping[1]) == 1
  assert len(mapping[0][0]) == 2
  assert len(mapping[1][0]) == 2
  assert abs(mapping[0][0]['Si'] - cos(2e0*pi*0e0/2.0)) < 1e-8
  assert abs(mapping[0][0]['Ge'] - cos(2e0*pi*1e0/2.0)) < 1e-8
  assert abs(mapping[1][0]['Si'] - cos(2e0*pi*0e0/2.0)) < 1e-8
  assert abs(mapping[1][0]['Ge'] - cos(2e0*pi*1e0/2.0)) < 1e-8

  lattice = binary.zinc_blende()
  lattice[1].type = ['Si', 'Ge', 'C']
  a = Cluster(lattice)
  mapping = a.occupation_mapping()
  assert len(mapping) == len(lattice)
  assert mapping[0] is None
  assert len(mapping[1]) == len(lattice[1].type) - 1
  assert len(mapping[1][0]) == 3
  assert abs(mapping[1][0]['C']  - cos(2e0*pi*0e0/3.0))
  assert abs(mapping[1][0]['Si'] - cos(2e0*pi*1e0/3.0))
  assert abs(mapping[1][0]['Ge'] - cos(2e0*pi*2e0/3.0))
  assert len(mapping[1][1]) == 3
  assert abs(mapping[1][0]['C']  - sin(2e0*pi*0e0/3.0))
  assert abs(mapping[1][0]['Si'] - sin(2e0*pi*1e0/3.0))
  assert abs(mapping[1][0]['Ge'] - sin(2e0*pi*2e0/3.0))

def test_onsite():
  from numpy import dot
  from random import choice
  from lada.crystal import binary, supercell
  from lada.ce import Cluster

  lattice = binary.zinc_blende()
  for atom in lattice: atom.type = ['Si', 'Ge']

  structure = binary.zinc_blende()
  for atom in structure: atom.type = 'Si'

  a = Cluster(lattice)

  # Empty cluster first
  assert abs(a(structure) - len(structure)) < 1e-8
  for i in xrange(10):
    superstructure = supercell(lattice, dot(lattice.cell, get_cell()))
    for atom in superstructure: atom.type = choice(atom.type)
    assert abs(a(superstructure) - len(superstructure)) < 1e-8

  # Try on-site cluster.
  # loop over random supercells.
  # PI should be number of proportional to number of each atomic type on each
  # site, or thereabouts
  mapping = a.occupation_mapping()
  for i in xrange(10):
    # create random superstructure
    superstructure = supercell(lattice, dot(lattice.cell, get_cell()))
    for atom in superstructure: atom.type = choice(atom.type)

    # now first and second site clusters
    for i, site in enumerate(lattice):
      # loop over flavors.
      types = [u.type for u in superstructure]
      for flavor in xrange(len(site.type) -1):
        a.spins = None
        a.add_spin(site.pos, flavor)
        s = 0e0
        for t in site.type:
          s += float(types.count(t)) * mapping[i][flavor][t]
        assert abs(a(superstructure) - s) < 1e-8

def test_symmetrized():
  from numpy import all, dot, any
  from numpy.random import randint
  from random import choice
  from lada.ce import Cluster
  from lada.crystal import binary, neighbors

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']

  a = Cluster(lattice)
  a.add_spin(lattice[0].pos, 0)
  a.add_spin(lattice[0].pos + [0.0, -0.5, -0.5], 0)
  a._create_symmetrized()

  assert len(a._symmetrized) == 2 * 12
  for i, (atom, vec, d) in enumerate(neighbors(lattice, 16, lattice[0].pos)):
    if i < 4: continue
    b = Cluster(lattice)
    b.add_spin(lattice[0].pos, 0)
    b.add_spin(vec, 0)
    assert any(all(b.spins == u) for u in a._symmetrized)
  for i, (atom, vec, d) in enumerate(neighbors(lattice, 16, lattice[1].pos)):
    if i < 4: continue
    b = Cluster(lattice)
    b.add_spin(lattice[1].pos, 0)
    b.add_spin(vec, 0)
    assert any(all(b.spins == u) for u in a._symmetrized)

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = Cluster(lattice)
  a.add_spin(lattice[1].pos, 1)
  a.add_spin(lattice[1].pos + [1.0, 0, 0], 0)
  a._create_symmetrized()

  assert len(a._symmetrized) ==  6
  for i, (atom, vec, d) in enumerate(neighbors(lattice, 24, lattice[0].pos)):
    if i < 16: continue
    b = Cluster(lattice)
    b.add_spin(lattice[1].pos, 1)
    b.add_spin(vec, 0)
    assert any(all(b.spins == u) for u in a._symmetrized)

def test_random():
  from numpy import dot
  from numpy.random import randint
  from random import choice
  from lada.ce import Cluster
  from lada.crystal import binary, supercell

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  # now try random clusters with cell and their supercells
  for i in xrange(10):
    # random cluster
    a = Cluster(lattice)
    site = choice([0, 1])
    flavor = 0
    if site == 1: flavor = choice([0, 1])
    a.add_spin(lattice[site].pos, flavor)
    for j in xrange(4): 
      site = choice([0, 1])
      flavor = 0
      if site == 1: flavor = choice([0, 1])
      pos = lattice[site].pos + dot(lattice.cell, randint(4, size=(3,))-2)
      try: a.add_spin(pos, flavor)
      except ValueError: pass

    # random structure
    structure = supercell(lattice, dot(lattice.cell, get_cell(5)))
    for atom in structure: atom.type = choice(lattice[atom.site].type)

    # random supercell
    for j in xrange(5):
      sp = supercell(structure, dot(structure.cell, get_cell(5)))
      for atom in sp: atom.site = structure[atom.site].site
      assert abs(a(sp) - len(sp) / len(structure) * a(structure)) < 1e-8

  

if __name__ == '__main__': 
  test_occmap()
  test_addspins()
  test_cmp()
  test_spins_are_sorted()
  test_onsite()
  test_symmetrized()
  test_random()
