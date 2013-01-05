def test_magmom():
  from pickle import loads, dumps
  from pylada.crystal.cppwrappers import Structure
  from pylada.vasp import Vasp

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

  vasp = Vasp()
  vasp.magmom = None
  assert vasp.magmom is None
  assert vasp._input['magmom'].keyword == 'MAGMOM'
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure) is None

  # ispin == 1
  vasp.magmom = True
  vasp.ispin = 1
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure) is None
  # ispins == 2, magmom == False
  vasp.ispin = 2
  vasp.magmom = None
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure) is None
  vasp.magmom = False
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure) is None
  # now for real print
  vasp.magmom = True
  assert 'MAGMOM' in vasp._input['magmom'].output_map(vasp=vasp, structure=structure)
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure)['MAGMOM'] == '4*-1.00 2*0.00 8*0.00'
  # now print a string directly.
  vasp.magmom = 'hello'
  assert vasp.magmom == 'hello'
  assert 'MAGMOM' in vasp._input['magmom'].output_map(vasp=vasp, structure=structure)
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure)['MAGMOM'] == 'hello'

  # check repr
  assert repr(vasp._input['magmom']) == "Magmom(value='hello')"
  # check pickling
  assert repr(loads(dumps(vasp._input['magmom']))) == "Magmom(value='hello')"

  # more tests.
  for atom, mag in zip(structure, [1, -1, 1, 1]):
    if atom.type == 'Mg': atom.magmom = mag
  for atom, mag in zip(structure, [0.5, 0.5]):
    if atom.type == 'Al': atom.magmom = mag

  vasp.magmom = True
  assert 'MAGMOM' in vasp._input['magmom'].output_map(vasp=vasp, structure=structure)
  assert vasp._input['magmom'].output_map(vasp=vasp, structure=structure)['MAGMOM'] == '1.00 -1.00 2*1.00 2*0.00 8*0.00'

if __name__ == "__main__":
  test_magmom()
