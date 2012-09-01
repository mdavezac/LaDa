def test_J0():
  from lada.ce import cluster_factory
  from lada.crystal import binary

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, J0=True)
  assert len(a) == 1
  assert a[0].order == 0

def test_J1():
  from numpy import all
  from lada.ce import cluster_factory
  from lada.crystal import binary

  # test multi-lattice with different occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, J1=True)
  assert len(a) == 2
  assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
  assert a[0].spins['sublattice'] == 0
  assert a[1].spins['sublattice'] == 1

  # test multi-lattice with same occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge']

  a = cluster_factory(lattice, J1=True)
  assert len(a) == 1
  assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
  assert a[0].spins['sublattice'] == 0

def test_B2():
  from numpy import all
  from lada.ce import cluster_factory
  from lada.crystal import binary

  # test multi-lattice with different occupations. 
  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  a = cluster_factory(lattice, B2=1)
  for u in a: print u
# assert len(a) == 2
# assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
# assert a[0].spins['sublattice'] == 0
# assert a[1].spins['sublattice'] == 1

if __name__ == '__main__': 
  test_J0()
  test_J1()
  test_B2()
