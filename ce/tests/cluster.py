def test():
  from lada.crystal import binary
  from lada.ce import Cluster

  lattice = binary.zinc_blende()
  for atom in lattice: atom.type = ['Si', 'Ge']

  structure = binary.zinc_blende()
  for atom in structure: atom.type = 'Si'

  a = Cluster(lattice)
  print a(structure)


if __name__ == '__main__': test()
