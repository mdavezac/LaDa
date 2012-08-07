""" Test for supercell function. """
def test_supercell():
  """ Simple supercell test. """
  from numpy import identity, abs, all, dot
  from numpy.linalg import inv
  from lada.crystal.cppwrappers import supercell, Structure, are_periodic_images as api
  lattice = Structure( 0.0, 0.5, 0.5,
                       0.5, 0.0, 0.5,
                       0.5, 0.5, 0.0, scale=2.0, m=True ) \
                     .add_atom(0, 0, 0, "As")           \
                     .add_atom(0.25, 0.25, 0.25, ['In', 'Ga'], m = True)
  result = supercell(lattice, dot(lattice.cell, [ [-1, 1, 1],
                                [1, -1, 1], 
                                [1, 1, -1] ] ) )
  assert all(abs(result.cell - identity(3)) < 1e-8)
  assert abs(result.scale - 2) < 1e-8
  assert getattr(result, 'm', False) 
  assert all(abs(result[0].pos - [0.00, 0.00, 0.00]) < 1e-8) and result[0].type == "As" \
         and getattr(result[0], 'site', -1) == 0 and api(result[0].pos, lattice[0].pos, inv(lattice.cell))
  assert all(abs(result[1].pos - [0.25, 0.25, 0.25]) < 1e-8) and result[1].type == ["In", "Ga"] \
         and getattr(result[1], 'm', False) and getattr(result[1], 'site', -1) == 1 \
         and api(result[1].pos, lattice[1].pos, inv(lattice.cell))
  assert all(abs(result[2].pos - [0.50, 0.00, 0.50]) < 1e-8) and result[2].type == "As" \
         and getattr(result[2], 'site', -1) == 0 and api(result[2].pos, lattice[0].pos, inv(lattice.cell))
  assert all(abs(result[3].pos - [0.75, 0.25, 0.75]) < 1e-8) and result[3].type == ["In", "Ga"] \
         and getattr(result[3], 'm', False) and getattr(result[3], 'site', -1) == 1 \
         and api(result[3].pos, lattice[1].pos, inv(lattice.cell))
  assert all(abs(result[4].pos - [0.50, 0.50, 0.00]) < 1e-8) and result[4].type == "As" \
         and getattr(result[4], 'site', -1) == 0 and api(result[4].pos, lattice[0].pos, inv(lattice.cell))
  assert all(abs(result[5].pos - [0.75, 0.75, 0.25]) < 1e-8) and result[5].type == ["In", "Ga"] \
         and getattr(result[5], 'm', False) and getattr(result[5], 'site', -1) == 1 \
         and api(result[5].pos, lattice[1].pos, inv(lattice.cell))
  assert all(abs(result[6].pos - [0.00, 0.50, 0.50]) < 1e-8) and result[6].type == "As" \
         and getattr(result[6], 'site', -1) == 0 and api(result[6].pos, lattice[0].pos, inv(lattice.cell))
  assert all(abs(result[7].pos - [0.25, 0.75, 0.75]) < 1e-8) and result[7].type == ["In", "Ga"] \
         and getattr(result[7], 'm', False) and getattr(result[7], 'site', -1) == 1 \
         and api(result[7].pos, lattice[1].pos, inv(lattice.cell))
         
if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])

  test_supercell()
