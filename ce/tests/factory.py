def test():
  from numpy import all
  from lada.ce import cluster_factory
  from lada.ce.cluster import spin
  from lada.crystal import binary

  lattice = binary.zinc_blende()
  lattice[0].type = ['Si', 'Ge']
  lattice[1].type = ['Si', 'Ge', 'C']

  # Test J0
  a = cluster_factory(lattice, J0=True)
  assert len(a) == 1
  assert a[0].order == 0
  # Test J1
  a = cluster_factory(lattice, J1=True)
  assert len(a) == 3
  assert all(all(abs(u.spins['position']) < 1e-8) for u in a)
  assert a[0].spins['sublattice'] == 0
  assert a[0].spins['flavor'] == 0
  assert a[1].spins['sublattice'] == 1
  assert a[1].spins['flavor'] == 0
  assert a[2].spins['sublattice'] == 1
  assert a[2].spins['flavor'] == 1
  # Test J2
  a = cluster_factory(lattice, B2=2)
  for u in a: print u
  return
# assert len(a) == 6
# assert all(all(u.spins[0] == spin([0,0,0])) for u in a[:3])
# assert all(all(u.spins[0] == spin([0,0,0], 1)) for u in a[3:])
# assert all(a[0].spins[1] == spin([0,0,0], 1))
# assert all(a[1].spins[1] == spin([0,0,0], 1, 1))
# assert all(a[2].spins[1] == spin([-0.5,0.5,0]))
# assert all(a[3].spins[1] == spin([0.5,0, 0.5]))
# assert all(a[4].spins[1] == spin([-0.5, 0.5, 0], 1))
# assert all(a[5].spins[1] == spin([-0.5, 0.5, 0], 1, 1))
# # Test J3
  a = cluster_factory(lattice, B4=2)
  for u in a: print u
def test_flatten():
  from lada.ce.factory import _flatten
  assert _flatten([0,1], 3) == [0, 1, 3]
  assert _flatten([0,1], [2], 3, [[4,5], 6]) == range(7)
  

if __name__ == '__main__': 
  test_flatten()
  test()
