def test():
  from collections import namedtuple
  from pickle import loads, dumps
  from numpy import all, abs, array
  from lada.crystal import Structure
  from lada.vasp import Vasp
  from lada.vasp.specie import U, nlep
  import lada

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
  a = Vasp()
  Specie = namedtuple('Specie', ['U'])
  a.species = {'A': Specie([]), 'B': Specie([]), 'X': Specie([])}
  lada.vasp_has_nlep  = False

  o = a._input['ldau']
  d = {'LDAU': o.__class__}
  assert a.ldau == True
  assert o.output_map(vasp=a, structure=structure) is None
  assert eval(repr(o), d)._value == True
  assert eval(repr(o), d).keyword == 'LDAU'
  assert loads(dumps(o)).keyword == 'LDAU'
  assert loads(dumps(o))._value 

  # now disables U.
  a.species = {'A': Specie([U(2, 0, 0.5)]), 'B': Specie([]), 'X': Specie([])}
  a.ldau = False
  assert a.ldau == False
  assert o.output_map() is None
  # now prints U
  a.ldau = True
  map = o.output_map(vasp=a, structure=structure)
  assert map['LDAU'] == '.TRUE.'
  assert map['LDAUTYPE'] == '2'
  assert all(abs(array(map['LDUJ'].split(), dtype='float64')) < 1e-8)
  assert all(abs(array(map['LDUU'].split(), dtype='float64')-[0.5, 0, 0]) < 1e-8)
  assert all(abs(array(map['LDUL'].split(), dtype='float64')-[0, -1, -1]) < 1e-8)
  a.species = {'A': Specie([U(2, 0, 0.5)]), 'B': Specie([U(2, 1, 0.6)]), 'X': Specie([])}
  map = o.output_map(vasp=a, structure=structure)
  assert map['LDAU'] == '.TRUE.'
  assert map['LDAUTYPE'] == '2'
  assert all(abs(array(map['LDUJ'].split(), dtype='float64')) < 1e-8)
  assert all(abs(array(map['LDUU'].split(), dtype='float64')-[0.5, 0, 0.6]) < 1e-8)
  assert all(abs(array(map['LDUL'].split(), dtype='float64')-[0, -1, 1]) < 1e-8)
  

  # now tries NLEP
  lada.vasp_has_nlep = True
  a.species = {'A': Specie([U(2, 0, 0.5)]), 'B': Specie([U(2, 0, -0.5), nlep(2, 1, -1.0)]), 'X': Specie([])}
  a.ldau = False
  assert a.ldau == False
  assert o.output_map() is None
  a.ldau = True
  map = o.output_map(vasp=a, structure=structure)
  assert map['LDAU'] == '.TRUE.'
  assert map['LDAUTYPE'] == '2'
  assert all(abs(array(map['LDUL1'].split(), dtype='float64')-[0, -1, 0]) < 1e-8)
  assert all(abs(array(map['LDUU1'].split(), dtype='float64')-[0.5, 0, -0.5]) < 1e-8)
  assert all(abs(array(map['LDUJ1'].split(), dtype='float64')-[0, 0, 0]) < 1e-8)
  assert all(abs(array(map['LDUO1'].split(), dtype='float64')-[1, 1, 1]) < 1e-8)
  assert all(abs(array(map['LDUL2'].split(), dtype='float64')-[-1, -1, 1]) < 1e-8)
  assert all(abs(array(map['LDUU2'].split(), dtype='float64')-[0, 0, -1.0]) < 1e-8)
  assert all(abs(array(map['LDUJ2'].split(), dtype='float64')-[0, 0, 0]) < 1e-8)
  assert all(abs(array(map['LDUO2'].split(), dtype='float64')-[1, 1, 2]) < 1e-8)

  a.species = {'A': Specie([U(2, 0, 0.5)]), 'B': Specie([U(2, 0, -0.5), nlep(2, 2, -1.0, -3.0)]), 'X': Specie([])}
  a.ldau = True
  map = o.output_map(vasp=a, structure=structure)
  assert map['LDAU'] == '.TRUE.'
  assert map['LDAUTYPE'] == '2'
  assert all(abs(array(map['LDUL1'].split(), dtype='float64')-[0, -1, 0]) < 1e-8)
  assert all(abs(array(map['LDUU1'].split(), dtype='float64')-[0.5, 0, -0.5]) < 1e-8)
  assert all(abs(array(map['LDUJ1'].split(), dtype='float64')-[0, 0, 0]) < 1e-8)
  assert all(abs(array(map['LDUO1'].split(), dtype='float64')-[1, 1, 1]) < 1e-8)
  assert all(abs(array(map['LDUL2'].split(), dtype='float64')-[-1, -1, 2]) < 1e-8)
  assert all(abs(array(map['LDUU2'].split(), dtype='float64')-[0, 0, -1.0]) < 1e-8)
  assert all(abs(array(map['LDUJ2'].split(), dtype='float64')-[0, 0, -3.0]) < 1e-8)
  assert all(abs(array(map['LDUO2'].split(), dtype='float64')-[1, 1, 3]) < 1e-8)


if __name__ == "__main__":
  test()

