def test_iterator():
  """ Test defect iterator. """
  from numpy import ones, logical_not
  from lada.enum.defects import Iterator

  size = 8
  mask = ones(size, dtype='bool')
  template = ones(size, dtype='int16')
  mask0 = mask.copy()
  mask0[::4] = False
  mask0[1::4] = False
  mask2 = mask.copy()
  mask2[2:] = False


  results = [[4, 3, 2, 2, 3, 1, 1, 1], [3, 4, 2, 2, 3, 1, 1, 1], [4, 3, 2, 1,
    3, 1, 2, 1], [3, 4, 2, 1, 3, 1, 2, 1], [4, 3, 2, 1, 3, 1, 1, 2], [3, 4, 2,
      1, 3, 1, 1, 2], [4, 3, 1, 2, 3, 1, 2, 1], [3, 4, 1, 2, 3, 1, 2, 1], [4,
        3, 1, 2, 3, 1, 1, 2], [3, 4, 1, 2, 3, 1, 1, 2], [4, 3, 1, 1, 3, 1, 2,
          2], [3, 4, 1, 1, 3, 1, 2, 2], [4, 3, 2, 2, 1, 3, 1, 1], [3, 4, 2, 2,
            1, 3, 1, 1], [4, 3, 2, 1, 1, 3, 2, 1], [3, 4, 2, 1, 1, 3, 2, 1],
          [4, 3, 2, 1, 1, 3, 1, 2], [3, 4, 2, 1, 1, 3, 1, 2], [4, 3, 1, 2, 1,
            3, 2, 1], [3, 4, 1, 2, 1, 3, 2, 1], [4, 3, 1, 2, 1, 3, 1, 2], [3,
              4, 1, 2, 1, 3, 1, 2], [4, 3, 1, 1, 1, 3, 2, 2], [3, 4, 1, 1, 1,
                3, 2, 2], [4, 1, 2, 2, 3, 3, 1, 1], [1, 4, 2, 2, 3, 3, 1, 1],
              [4, 1, 2, 1, 3, 3, 2, 1], [1, 4, 2, 1, 3, 3, 2, 1], [4, 1, 2, 1,
                3, 3, 1, 2], [1, 4, 2, 1, 3, 3, 1, 2], [4, 1, 1, 2, 3, 3, 2,
                  1], [1, 4, 1, 2, 3, 3, 2, 1], [4, 1, 1, 2, 3, 3, 1, 2], [1,
                    4, 1, 2, 3, 3, 1, 2], [4, 1, 1, 1, 3, 3, 2, 2], [1, 4, 1,
                      1, 3, 3, 2, 2]]

  iterator = Iterator(template, (1, 4, mask2), (2, 2, mask0), (2, 3, logical_not(mask0)))
  for i, x in enumerate(iterator):
    assert all(x == results[i])
# def test_defects():
#   from lada.crystal.bravais import fcc
#   from lada.enum.defects import defects
#   lattice = fcc()
#   lattice[0].type = 'Zr', 'Ti'
#   lattice.add_atom(0.25, 0.25, 0.25, 'O', 'V')
#   lattice.add_atom(0.75, 0.75, 0.75, 'O', 'V')


#   defects(lattice, 8, {'V': 2, 'Ti': 2})


if __name__ == '__main__':
  test_iterator()
# test_defects()
