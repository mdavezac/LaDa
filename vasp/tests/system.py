def test_system():
  from lada.crystal.cppwrappers import Structure
  from lada.vasp import Vasp
  a = Vasp()
  b = Structure()

  assert a.system is None
  assert a._input['system'].keyword == 'system'
  assert a._input['system'].output_map(vasp=a, structure=b) is None

  a.system = 'system'
  assert a.system == 'system'
  assert 'system' in a._input['system'].output_map(vasp=a, structure=b)
  assert a._input['system'].output_map(vasp=a, structure=b)['system'] == 'system'

  b.name = 'hello'
  assert 'system' in a._input['system'].output_map(vasp=a, structure=b)
  assert a._input['system'].output_map(vasp=a, structure=b)['system'] == 'system: hello'

  a.system = None
  assert a.system is None
  assert 'system' in a._input['system'].output_map(vasp=a, structure=b)
  assert a._input['system'].output_map(vasp=a, structure=b)['system'] == 'hello'

  a.system = None
  assert a.system is None
  assert 'system' in a._input['system'].output_map(vasp=a, structure=b)
  assert a._input['system'].output_map(vasp=a, structure=b)['system'] == 'hello'

if __name__ == "__main__": 
  test_system()
