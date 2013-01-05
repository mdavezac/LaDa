def test(path):
  from pickle import loads, dumps
  from pylada.vasp import Vasp
  from pylada.crystal import Structure

  structure = Structure([[0, 0.5, 0.5],[0.5, 0, 0.5], [0.5, 0.5, 0]], scale=5.43, name='has a name')\
                       .add_atom(0,0,0, "Si")\
                       .add_atom(0.25,0.25,0.25, "Si")
  a = Vasp()
  a.add_specie = "Si", "{0}/pseudos/Si".format(path)
  assert a.extraelectron is None
  assert a._input['extraelectron'].output_map() is None
  assert a._input['nelect'].output_map() is None
  a.extraelectron = 0
  assert a.extraelectron == 0
  assert a.nelect is None
  assert a._input['extraelectron'].output_map() is None
  assert a._input['nelect'].output_map() is None
  a.extraelectron = 1
  assert a.extraelectron == 1
  assert a.nelect is None
  assert 'nelect' in a._input['extraelectron'].output_map(vasp=a, structure=structure) 
  assert abs(float(a._input['extraelectron'].output_map(vasp=a, structure=structure)['nelect']) - 9.0) < 1e-8
  assert a._input['nelect'].output_map() is None
  a.nelect = 1
  a.extraelectron = -1
  assert a.extraelectron == -1
  assert a.nelect is None
  assert 'nelect' in a._input['extraelectron'].output_map(vasp=a, structure=structure) 
  assert abs(float(a._input['extraelectron'].output_map(vasp=a, structure=structure)['nelect']) - 7.0) < 1e-8
  assert a._input['nelect'].output_map() is None
  o = a._input['extraelectron']
  d = {'ExtraElectron': o.__class__}
  assert repr(eval(repr(o), d)) == repr(o)
  assert abs(float(eval(repr(o), d).output_map(vasp=a, structure=structure)['nelect']) - 7.0) < 1e-8
  assert repr(loads(dumps(o))) == repr(o)

  a.nelect = 8
  assert a.nelect == 8
  assert a.extraelectron is None
  assert 'nelect' in a._input['nelect'].output_map()
  assert abs(float(a._input['nelect'].output_map()['nelect']) - 8.0) < 1e-8
  assert a._input['extraelectron'].output_map() is None
  o = a._input['nelect']
  d = {'NElect': o.__class__}
  assert repr(eval(repr(o), d)) == repr(o)
  assert abs(float(eval(repr(o), d).output_map()['nelect']) - 8.0) < 1e-8
  assert repr(loads(dumps(o))) == repr(o)

if __name__ == "__main__":
  from sys import argv
  test(argv[1])
  
