def test_magmom():
  from collections import namedtuple
  from pickle import loads, dumps
  from lada.crystal.cppwrappers import Structure
  from lada.vasp.incar._params import Magmom

  u = 0.25
  x, y = u, 0.25-u
  structure = Structure([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]) \
                       .add_atom(5.000000e-01, 5.000000e-01, 5.000000e-01, "Mg") \
                       .add_atom(5.000000e-01, 2.500000e-01, 2.500000e-01, "Mg") \
                       .add_atom(2.500000e-01, 5.000000e-01, 2.500000e-01, "Mg") \
                       .add_atom(2.500000e-01, 2.500000e-01, 5.000000e-01, "Mg") \
                       .add_atom(8.750000e-01, 8.750000e-01, 8.750000e-01, "Al") \
                       .add_atom(1.250000e-01, 1.250000e-01, 1.250000e-01, "Al") \
                       .add_atom(     x,     x,     x, "O") \
                       .add_atom(     x,     y,     y, "O") \
                       .add_atom(     y,     x,     y, "O") \
                       .add_atom(     y,     y,     x, "O") \
                       .add_atom(    -x,    -x,    -x, "O") \
                       .add_atom(    -x,    -y,    -y, "O") \
                       .add_atom(    -y,    -x,    -y, "O") \
                       .add_atom(    -y,    -y,    -x, "O") 
  
  for atom in structure:
    if atom.type == 'Mg': atom.magmom = -1e0

  Vasp = namedtuple('Vasp', ['ispin'])
  vasp = Vasp(1)

  # ispin == 1
  magmom = Magmom(True)
  assert magmom.incar_string(vasp=vasp, structure=structure) is None
  # ispins == 2, magmom == False
  magmom.value = False
  vasp = Vasp(2)
  assert magmom.incar_string(vasp=vasp, structure=structure) is None
  # now for real print
  magmom.value = True
  assert magmom.incar_string(vasp=vasp, structure=structure) == 'MAGMOM = 4*-1.00 2*0.00 8*0.00'
  # now print a string directly.
  magmom.value = 'hello'
  assert magmom.incar_string(vasp=vasp, structure=structure) == 'MAGMOM = hello'

  # check repr
  assert repr(magmom) == "Magmom('hello')"
  # check pickling
  assert repr(loads(dumps(magmom))) == "Magmom('hello')"

  # more tests.
  magmom.value = True
  for atom, mag in zip(structure, [1, -1, 1, 1]):
    if atom.type == 'Mg': atom.magmom = mag
  for atom, mag in zip(structure, [0.5, 0.5]):
    if atom.type == 'Al': atom.magmom = mag

  assert magmom.incar_string(vasp=vasp, structure=structure) \
           == 'MAGMOM = 1.00 -1.00 2*1.00 2*0.00 8*0.00'



if __name__ == "__main__":
  from sys import argv, path 
  from numpy import array
  if len(argv) > 0: path.extend(argv[1:])
  
  test_magmom()

