def test():
  from collections import namedtuple
  from pickle import loads, dumps
  from pylada.crystal.cppwrappers import Structure
  from pylada.vasp.incar._params import ExtraElectron

  u = 0.25
  x, y = u, 0.25-u
  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
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
  
  for atom in structure:
    if atom.type == 'Mg': atom.magmom = -1e0

  Vasp = namedtuple('Vasp', ['species'])
  Specie = namedtuple('Specie', ['valence'])
  vasp = Vasp({'A': Specie(5), 'B': Specie(10), 'X': Specie(2)})

  assert 4*5+2*10+2*8 == ExtraElectron(0).nelectrons(vasp, structure)
  assert ExtraElectron(0).incar_string(vasp=vasp, structure=structure)\
           == "# NELECT = 56.0 Charge neutral system"
  assert ExtraElectron(1).incar_string(vasp=vasp, structure=structure)\
           == "NELECT = 57.0  # negatively charged system (-1)"
  assert ExtraElectron(-1).incar_string(vasp=vasp, structure=structure)\
           == "NELECT = 55.0  # positively charged system (+1)"
  assert repr(ExtraElectron()) == 'ExtraElectron(0)'
  assert repr(loads(dumps(ExtraElectron()))) == 'ExtraElectron(0)'
  assert repr(ExtraElectron(-1)) == 'ExtraElectron(-1)'
  assert repr(loads(dumps(ExtraElectron(-1)))) == 'ExtraElectron(-1)'


if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test()

