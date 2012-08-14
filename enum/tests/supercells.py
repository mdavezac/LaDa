def test():
  from lada.crystal.bravais import fcc
  from lada.enum import supercells

  # First, a lattice without symmetry
  lattice = fcc()
  lattice.cell = [[1.0, 0.2, 0], [0, 0.9, -0.1], [0, 0, 0.8]]
  scs = supercells(lattice, range(17))
  results = [1, 7, 13, 35, 31, 91, 57, 155, 130, 217, 133, 455, 183, 399, 403, 651]
  for r, s in zip(results, scs.itervalues()):
    assert r == len(s)
  # Then, an fcc lattice.
  lattice = fcc()
  scs = supercells(lattice, range(11))
  results = [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
  for r, s in zip(results, scs.itervalues()):
    assert r, len(s)
  

if __name__ == '__main__':
  test()
