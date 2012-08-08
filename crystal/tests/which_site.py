""" Checks that map_sites is correct. """
def test_b5(u):
  """ Test b5 space-group and equivalents """
  from random import randint, random
  from numpy import array
  from numpy.linalg import det
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
  for i, atom in enumerate(lattice):
    print which_site(atom.pos, lattice)
    assert which_site(atom.pos, lattice) == i

if __name__ == "__main__":
  from sys import argv, path 
  if len(argv) > 0: path.extend(argv[1:])
  
  test_b5(0.25)
  test_b5(0.36)
