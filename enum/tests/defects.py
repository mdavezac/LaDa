def test_iterator():
  """ Test defect iterator. """
  from numpy import ones, logical_not
  from lada.enum.defects import Iterator

  size = 8
  mask = ones(size, dtype='bool')
  mask0 = mask.copy()
  mask0[::4] = False
  mask0[1::4] = False
  mask2 = mask.copy()
  mask2[2:] = False


  results = [[4, 3, 2, 2, 3, 0, 0, 0], [3, 4, 2, 2, 3, 0, 0, 0], [4, 3, 2, 0,
    3, 0, 2, 0], [3, 4, 2, 0, 3, 0, 2, 0], [4, 3, 2, 0, 3, 0, 0, 2], [3, 4, 2,
      0, 3, 0, 0, 2], [4, 3, 0, 2, 3, 0, 2, 0], [3, 4, 0, 2, 3, 0, 2, 0], [4,
        3, 0, 2, 3, 0, 0, 2], [3, 4, 0, 2, 3, 0, 0, 2], [4, 3, 0, 0, 3, 0, 2,
          2], [3, 4, 0, 0, 3, 0, 2, 2], [4, 3, 2, 2, 0, 3, 0, 0], [3, 4, 2, 2,
            0, 3, 0, 0], [4, 3, 2, 0, 0, 3, 2, 0], [3, 4, 2, 0, 0, 3, 2, 0],
          [4, 3, 2, 0, 0, 3, 0, 2], [3, 4, 2, 0, 0, 3, 0, 2], [4, 3, 0, 2, 0,
            3, 2, 0], [3, 4, 0, 2, 0, 3, 2, 0], [4, 3, 0, 2, 0, 3, 0, 2], [3,
              4, 0, 2, 0, 3, 0, 2], [4, 3, 0, 0, 0, 3, 2, 2], [3, 4, 0, 0, 0,
                3, 2, 2], [4, 0, 2, 2, 3, 3, 0, 0], [0, 4, 2, 2, 3, 3, 0, 0],
              [4, 0, 2, 0, 3, 3, 2, 0], [0, 4, 2, 0, 3, 3, 2, 0], [4, 0, 2, 0,
                3, 3, 0, 2], [0, 4, 2, 0, 3, 3, 0, 2], [4, 0, 0, 2, 3, 3, 2,
                  0], [0, 4, 0, 2, 3, 3, 2, 0], [4, 0, 0, 2, 3, 3, 0, 2], [0,
                    4, 0, 2, 3, 3, 0, 2], [4, 0, 0, 0, 3, 3, 2, 2], [0, 4, 0,
                      0, 3, 3, 2, 2]]

  iterator = Iterator(size, (1, 4, mask2), (2, 2, mask0), (2, 3, logical_not(mask0)))
  for i, x in enumerate(iterator):
    assert all(x == results[i])
  assert i == len(results) - 1
  iterator.reset()
  doiter = False
  for i, x in enumerate(iterator):
    doiter = True
    assert all(x == results[i])
  assert doiter
  assert i == len(results) - 1

def test_defects():
  from lada.crystal.bravais import fcc
  from lada.enum.defects import defects
  lattice = fcc()
  lattice[0].type = 'Zr', 'Ti'
  lattice.add_atom(0.25, 0.25, 0.25, 'O', 'A')
  lattice.add_atom(0.75, 0.75, 0.75, 'O', 'A')


  for i, (x, hft, hermite) in enumerate(defects(lattice, 32, {'A': 2, 'Ti': 2})):
    continue
  print i


if __name__ == '__main__':
  test_iterator()
  test_defects()
