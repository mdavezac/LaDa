def test_tree():
  from numpy import all, array, dot, sum, any
  from lada.crystal import binary, supercell
  from lada.vff import Functional
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = Functional.build_tree(a, overlap=0.5)
  
  for center in b:
    positions = []
    for i, (bond, vector) in enumerate(center):
      position = bond.pos + dot(a.cell, vector)
      assert abs(sum((position - center.pos)**2) - 0.25*0.25*3) < 1e-8
      assert all( [any(abs(array(p) - position[None, :]) > 1e-8) for p in positions] )
      positions.append(position)
    assert i == 3

def test_disorder(lim=8):
  from numpy import all, array, dot, sum, any, identity
  from numpy.random import random, randint
  from numpy.linalg import det
  from lada.crystal import binary, supercell
  from lada.vff import Functional


  lattice = binary.zinc_blende()
  for i in xrange(10):
    cell = randint(-lim, lim, (3,3))
    while det(cell) == 0: cell = randint(-lim, lim, (3,3))
    a = supercell(lattice, dot(lattice.cell, cell))
  
    b = Functional.build_tree(a, overlap=0.8)
    ids = [id(node.center) for node in b]
    connections = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
  
    epsilon = random((3,3)) * 0.1
    epsilon = epsilon + epsilon.T
    a.cell += dot(epsilon, a.cell)
    for atom in a: atom.pos += dot(epsilon, atom.pos)
    
    b = Functional.build_tree(a, overlap=0.8)
    c = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
    assert all(connections == c)
    
    b = Functional.build_tree(a, overlap=0.8)
    for atom in a: atom.pos += random(3) * 0.05 - 0.025
    c = array([ sorted([ids.index(id(n.center)) for n, v in node])
                          for node in b ])
    assert all(connections == c)
  return a

def test_single_counting():
  from lada.crystal import binary, supercell
  from lada.vff import Functional
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[4, 0, 0], [0, 2, 0], [0, 0, 1]])
  b = Functional.build_tree(a, overlap=0.5)
  
  n = 0
  for center in b:
    positions = []
    for endpoint, vector in center.sc_bond_iter():
      n += 1
      for other, v in endpoint.sc_bond_iter(): assert other is not center
      assert id(center) in [id(c) for c, v in endpoint]
  assert n == 2 * len(a)

if __name__ == '__main__':
  test_tree()
  test_single_counting()
# structure = test_disorder()
