def test():
  from numpy import abs, all
  from pylada.dftcrystal.parse import parse
  from pylada.dftcrystal import Molecule
  string = '\nMOLECULE\n1\n2\n1 0.0 -2.91352558499E-15 12.3754696347\n'        \
           '1 0.0 -2.91352558499E-15 13.1454696347\nEND MOLECULE'
  tree = parse(string)['']
  molecule = Molecule()
  molecule.read_input(tree['MOLECULE'])
  assert molecule.symmgroup == 1
  assert len(molecule) == 0
  assert len(molecule.atoms) == 2
  assert molecule.atoms[0].type == 'H'
  assert molecule.atoms[1].type == 'H'
  assert all(abs(molecule.atoms[0].pos - [0, 0, 12.3754696347]) < 1e-6)
  assert all(abs(molecule.atoms[1].pos - [0, 0, 13.1454696347]) < 1e-6)

if __name__ == '__main__': test()
