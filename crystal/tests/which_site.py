""" Checks that map_sites is correct. """
def test_b5(u):
  """ Test b5 space-group and equivalents """
  from numpy import dot
  from numpy.random import randint
  from lada.crystal import Structure, which_site

  x, y = u, 0.25-u
  lattice = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
                     .add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "A") \
                     .add_atom(5.000000e-01, 2.500000e-01, 2.500000e-01, "A") \
                     .add_atom(2.500000e-01, 5.000000e-01, 2.500000e-01, "A") \
                     .add_atom(2.500000e-01, 2.500000e-01, 5.000000e-01, "A") \
                     .add_atom(8.750000e-01, 8.750000e-01, 8.750000e-01, "B") \
                     .add_atom(1.250000e-01, 1.250000e-01, 1.250000e-01, "B") \
                     .add_atom(     x,     x,     x, "X") \
                     .add_atom(     x,     y,     y, "X") \
                     .add_atom(     y,     x,     y, "X") \
                     .add_atom(     y,     y,     x, "X") \
                     .add_atom(    -x,    -x,    -x, "X") \
                     .add_atom(    -x,    -y,    -y, "X") \
                     .add_atom(    -y,    -x,    -y, "X") \
                     .add_atom(    -y,    -y,    -x, "X") 
  assert which_site(lattice[6].pos + [0.5, -0.5, 2], lattice) == 6
  for i, atom in enumerate(lattice):
    assert which_site(atom.pos, lattice) == i
    for j in xrange(10):
      newpos = dot(lattice.cell, randint(10, size=(3,)) - 5)
      assert which_site(atom.pos + newpos, lattice) == i, (atom.pos, newpos, i)

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  test_b5(0.25)
  test_b5(0.36)
