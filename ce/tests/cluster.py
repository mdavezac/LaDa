def get_cell(n=5):
  from numpy.random import randint
  from numpy.linalg import det
  cell = randint(2*n, size=(3,3)) - n
  while abs(det(cell)) < 1e-8:
    cell = randint(2*n, size=(3,3)) - n
  return cell

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

def test():
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
      types = [u.type for u in superstructure if u.site == i]
      for flavor in xrange(len(site.type) -1):
        a.spins = None
        a.add_spin(site.pos, flavor)
        s = 0e0
        for t in site.type:
          s += float(types.count(t)) * mapping[i][flavor][t]
        assert abs(a(superstructure) - s) < 1e-8
        print 'here'


if __name__ == '__main__': 
  test_occmap()
  test()
