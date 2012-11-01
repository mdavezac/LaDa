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

def test_disorder():
  from numpy import all, array, dot, sum, any, identity
  from numpy.random import random
  from lada.crystal import binary, supercell
  from lada.vff import Functional
  a = binary.zinc_blende()
  a = supercell(binary.zinc_blende(), [[2, 0, 0], [0, 2, 0], [0, 0, 1]])

  epsilon = random((3,3)) * 0.1
  epsilon = epsilon + epsilon.T
  a.cell += dot(epsilon, a.cell)
  for atom in a: atom.pos += dot(epsilon, atom.pos)
  
  b = Functional.build_tree(a, overlap=0.5)
  for center in b:
    positions = []
    for i, (bond, vector) in enumerate(center):
      position = bond.pos + dot(a.cell, vector)
      assert all( [any(abs(array(p) - position[None, :]) > 1e-8) for p in positions] )
      positions.append(position)
    print i
#   assert i == 3
  return a
#   print i
#   assert i == 3

  
# for center in b:
#   positions = []
#   for i, (bond, vector) in enumerate(center):
#     position = bond.pos + dot(a.cell, vector)
#     assert abs(sum((position - center.pos)**2) - 0.25*0.25*3) < 1e-8
#     assert all( [any(abs(array(p) - position[None, :]) > 1e-8) for p in positions] )
#     positions.append(position)
#   assert i == 3
if __name__ == '__main__':
  test_tree()
  structure = test_disorder()
