""" Checks neighbor routine. """
def check(structure, center, tolerance=1e-8):
  from numpy import abs, sqrt, all
  from lada.crystal.cppwrappers import neighbors
  
  # check we get the neighbors of zinc-blende.
  neighs = neighbors(structure, 46, center, tolerance);
  assert len(neighs) == 46
  for atom, diff, dist in neighs[:4]:
    assert abs(dist - sqrt(3.0) * 0.25) < tolerance
    assert all(abs(abs(diff) - 0.25) < tolerance)
  for atom, diff, dist in neighs[4: 16]:
    assert abs(dist - sqrt(2.0) * 0.5) < tolerance
    assert len([0 for u in diff if abs(u) < tolerance]) == 1
    assert len([0 for u in diff if abs(abs(u)-0.5) < tolerance]) == 2
  for atom, diff, dist in neighs[16: 28]:
    assert abs(dist - sqrt(0.75*0.75+2.*0.25*0.25)) < tolerance
    assert len([0 for u in diff if abs(abs(u)-0.25) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-0.75) < tolerance]) == 1
  for atom, diff, dist in neighs[28:34]:
    assert abs(dist - 1.) < tolerance
    assert len([0 for u in diff if abs(u) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-1.) < tolerance]) == 1
  for atom, diff, dist in neighs[34:]:
    assert abs(dist - sqrt(2*0.75*0.75+0.25*0.25)) < tolerance
    assert len([0 for u in diff if abs(abs(u)-0.75) < tolerance]) == 2
    assert len([0 for u in diff if abs(abs(u)-0.25) < tolerance]) == 1

  
  # check input using position rather than atom.
  neighs2 = neighbors(structure, 36, center.pos, tolerance);
  for (a, b, c), (d, e, f) in zip(neighs, neighs2):
    assert a is d
    assert all(abs(b-e)) < tolerance
    assert abs(c-f) < tolerance

  # check neighbor completeness.
  assert len(neighbors(structure, 2, center,tolerance)) == 4
  assert len(neighbors(structure, 4, center,tolerance)) == 4
  assert len(neighbors(structure, 6, center,tolerance)) == 16

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])

  from random import random
  from numpy import array
  from lada.crystal.cppwrappers import supercell, Structure

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
