def test():
  from lada.crystal.bravais import fcc
  from lada.enum.defects import defects
  lattice = fcc()
  lattice[0].type = 'Zr', 'Ti'
  lattice.add_atom(0.25, 0.25, 0.25, 'O', 'V')
  lattice.add_atom(0.75, 0.75, 0.75, 'O', 'V')


  defects(lattice, 8, {'V': 2, 'Ti': 2})


if __name__ == '__main__':
  test()
