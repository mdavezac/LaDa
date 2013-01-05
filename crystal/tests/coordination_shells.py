""" Checks coordination shell routine. """
def check(structure, center, tolerance=1e-8):
  from numpy import abs, sqrt, all
  from pylada.crystal.cppwrappers import coordination_shells
  
  # check we get the coordination_shells of zinc-blende.
  neighs = coordination_shells(structure, 5, center, tolerance);
  assert len(neighs) == 5
  assert len(neighs[0]) == 4
  for atom, diff, dist in neighs[0]:
    assert abs(dist - sqrt(3.0) * 0.25) < tolerance
    assert all(abs(abs(diff) - 0.25) < tolerance)
  assert len(neighs[1]) == 12
  for atom, diff, dist in neighs[1]:
    assert abs(dist - sqrt(2.0) * 0.5) < tolerance
    assert len([0 for u in diff if abs(u) < tolerance]) == 1
    assert len([0 for u in diff if abs(abs(u)-0.5) < tolerance]) == 2
  assert len(neighs[2]) == 12
  for atom, diff, dist in neighs[2]:
    assert abs(dist - sqrt(0.75*0.75+2.*0.25*0.25)) < tolerance
    assert len([0 for u in diff if abs(abs(u)-0.25) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-0.75) < tolerance]) == 1
  assert len(neighs[3]) == 6
  for atom, diff, dist in neighs[3]:
    assert abs(dist - 1.) < tolerance
    assert len([0 for u in diff if abs(u) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-1.) < tolerance]) == 1
  assert len(neighs[4]) == 12
  for atom, diff, dist in neighs[4]:
    assert abs(dist - sqrt(2*0.75*0.75+0.25*0.25)) < tolerance
    assert len([0 for u in diff if abs(abs(u)-0.75) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-0.25) < tolerance]) == 1

def check_against_neighbors(structure, tolerance=1e-8):
  from numpy import abs, sqrt, all
  from pylada.crystal.cppwrappers import coordination_shells, neighbors

  a = neighbors(structure, 150, [0,0,0], tolerance) 
  result = []
  fn = a[0][2]
  i = 0
  for atom, trans, dist in a:
    if abs(fn-dist) < tolerance: i+=1; continue
    result.append([i, fn])
    i = 1
    fn = dist

  b = coordination_shells(structure, 150, [0,0,0], tolerance) 
  for x, y in zip(result, b):
    assert len(y) == x[0]
    assert abs(y[0][2] - x[1]) < tolerance
  
if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])

  from random import random
  from numpy import array
  from pylada.crystal.cppwrappers import supercell, Structure

  lattice = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]]) \
                     .add_atom(0, 0, 0, "Si")                       \
                     .add_atom(0.25, 0.25, 0.25, "Ge")
  for atom in lattice: check(lattice, atom)

  structure = supercell(lattice, [1, 1, 0, -5, 2, 0, 0, 0, 1])
  for atom in structure: check(structure, atom)

  for atom in lattice: atom.pos += array([random(), random(), random()])*1e-4-5e-5 
  for atom in lattice: check(lattice, atom, 1e-2)

  for atom in structure: atom.pos += array([random(), random(), random()])*1e-4-5e-5 
  for atom in structure: check(structure, atom, 1e-2)

  check_against_neighbors(structure, 1e-2)

  lattice = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]]) \
                     .add_atom(0, 0, 0, "Si")
  check_against_neighbors(structure, 1e-8)
  lattice = Structure([[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]) \
                     .add_atom(0, 0, 0, "Si")
  check_against_neighbors(structure, 1e-8)
